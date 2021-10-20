require(reticulate)
require(imager)
if(!file.exists("MTprocess.py")){
	Stop("Python script missing")
}

get_pyids=function(allnets){
	#convert IDs into a python style


	allpyids=(lapply(allnets,function(x){
		ids=V(x)$name
		ids2=paste0("'",ids,"'")
		ids3=paste0(ids2,collapse=",")
		cmdstring=paste0("pyids=set({",ids3,"})")
		py_run_string(cmdstring)
		py$pyids
	}))
	return(allpyids)
}

get_com_pyids=function(coms){
	#convert community membership into python style

	allpycoms=lapply(coms,function(x){
		lapply(1:length(x),function(y){
			ids=x[[y]]
			ids2=paste0("'",ids,"'")
			ids3=paste0(ids2,collapse=",")
			cmdstring=paste0("pycoms=set({",ids3,"})")
			py_run_string(cmdstring)
			py$pycoms
		})
	})
	return(allpycoms)
}

do_track=function(allpycoms, allpyids,history=2){


	source_python("MTprocess.py")
	track = R_do_track(allpycoms, allpyids, history=history)
	return(track)
}

get_dc_membership=function(track){
	#get per timestep memberships of dynamic community from MajorTrack return
	dcmembership=lapply(track$individual_membership,function(x){
		unlist(x)
	})
	return(dcmembership)
}

add_dc_membership=function(allnets,dcmembership){
	#apply membership of dynamic community as node attribute
	allnets=lapply(1:length(allnets),function(x){
		V(allnets[[x]])$DC=dcmembership[[x]][match(V(allnets[[x]])$name,names(dcmembership[[x]]))]
		allnets[[x]]
	})
	return(allnets)
}

move_events_df=function(track,dcmembership=NULL,allremains=F){
	if(is.null(dcmembership)){
		dcmembership=get_dc_membership(track)
	}
	
	##Build dataframe of movement, splits, merges and remains between DCs.
	#get all split events
	newfromsplits=do.call(rbind,lapply(1:length(track$community_splits),function(slice){
		currslice=track$community_splits[[slice]]
		if(length(currslice)==0){
			return()
		}
		do.call(rbind,lapply(currslice,function(x){
			data.frame(slice=slice,parent=x[[1]],child=x[[2]],type=I("split"))
		}))	
	}))

	#get all merge events
	merges=do.call(rbind,lapply(1:length(track$community_merges),function(slice){
		currslice=track$community_merges[[slice]]
		if(length(currslice)==0){
			return()
		}
		do.call(rbind,lapply(currslice,function(x){
			data.frame(slice=slice,parent=x[[1]],child=x[[2]],type=I("merge"))
		}))
			
	}))
	

	#combine
	comorigins=rbind(newfromsplits,merges)

	#create unique move id
	comorigins$moveid=paste(comorigins$slice,comorigins$parent,comorigins$child)

	#reorder based on move ID, then slice
	comorigins=comorigins[order(comorigins$moveid),]
	comorigins=comorigins[order(comorigins$slice),]

	#where we have identical splits between and merges between groups set type to move
	comorigins$type[comorigins$moveid%in%comorigins$moveid[duplicated(comorigins$moveid)]]="move"
		
	#create remain category when parent and child are the same
	comorigins$type[comorigins$parent==comorigins$child]="remain"

	#remove the duplicate moves - also removes duplicate remains
	comorigins=comorigins[!duplicated(comorigins$moveid)|comorigins$type%in%c("split","merge"),]
	
	#get size of moves
	comorigins$size=sapply(1:nrow(comorigins),function(x){
		cevent=comorigins[x,]
		cslice=dcmembership[[cevent$slice]]
		#individuals in child event
		cmembers=names(cslice[cslice==cevent$child])
		pslice=dcmembership[[cevent$slice-1]]
		#individuals in parent event
		pmembers=names(pslice[pslice==cevent$parent])
		#individuals
		cmembers=cmembers[cmembers%in%pmembers] 
		length(cmembers)
	})
	
	if(allremains){
  	memdf=ind_membership_df(dcmembership)
  	newremains=lapply(unique(comorigins$slice),function(x){
  	  currinds=memdf[[1]][memdf[[1]]$timestep==x,]
  	  allgroups=unique(currinds$group)
  	  previnds=memdf[[1]][memdf[[1]]$timestep==x-1,]
  	  allprevgroups=unique(previnds$group)
  	  #get groups existing in previous timestep but not in comorigins
  	  missinggroups=allgroups[!allgroups%in%comorigins$child[comorigins$slice==x]&allgroups%in%allprevgroups]
  	  if(length(missinggroups)>0){
  	    data.frame(slice=x,parent=missinggroups,child=missinggroups,type="remain",moveid=paste(x,missinggroups,missinggroups),size=sapply(missinggroups,function(y){sum(currinds$group==y)}))
  	  }
  	  
   })
  	newremains=do.call(rbind,newremains)
  	comorigins=rbind(comorigins,newremains)
  	
  	#reorder based on move ID, then slice
  	comorigins=comorigins[order(comorigins$moveid),]
  	comorigins=comorigins[order(comorigins$slice),]
	}
	
	#ADD EMMIGRATION AND IMMIGRATION
	return(comorigins)
}

ind_membership_df=function(dcmembership){
	##A dataframe of individual DC membership per timestep
	memdf=do.call(rbind,lapply(1:length(dcmembership),function(x){
		data.frame(id=names(dcmembership[[x]]),timestep=x,group=dcmembership[[x]])
	}))
	
	##A matrix of individual group membership/presence over time
	allids=unique(memdf$id)
	allgroups=unique(memdf$group)
	alltimesteps=unique(memdf$timestep)
	memdf2=do.call(cbind,lapply(alltimesteps,function(x){
		cgroup=memdf$group[memdf$timestep==x][match(allids,memdf$id[memdf$timestep==x])]
		cgroup[!is.na(cgroup)]=cgroup[!is.na(cgroup)]
		cgroup[is.na(cgroup)]=NA
		return(cgroup)
	}))
	row.names(memdf2)=allids
	colnames(memdf2)=alltimesteps
	return(list(memdf1=memdf,memdf2=memdf2))
}

community_lifespans=function(track){
  unlist(track$community_lifespans)
}

get_DC_names=function(track,inputnames,timestep){
  track$comm_group_members[[timestep]]
  allnames=lapply(1:length(track$comm_group_members[[timestep]]),function(x){
    data.frame(DC=rep(names(track$comm_group_members[[timestep]])[[x]],length(track$comm_group_members[[timestep]][[x]])),
    cluster=track$comm_group_members[[timestep]][[x]])           
  })
  allnames=do.call(rbind,allnames)
  return(allnames$DC[match(inputnames,allnames$cluster)])     
}

get_cluster_names=function(track,inputnames,timestep){
  track$comm_group_members[[timestep]]
  allnames=lapply(1:length(track$comm_group_members[[timestep]]),function(x){
    data.frame(DC=rep(names(track$comm_group_members[[timestep]])[[x]],length(track$comm_group_members[[timestep]][[x]])),
               cluster=track$comm_group_members[[timestep]][[x]])           
  })
  allnames=do.call(rbind,allnames)
  return(allnames$cluster[match(inputnames,allnames$DC)])   
}

get_similarities=function(track){
  gs=track$group_similarities
  
  allsim=data.frame()
  for(i in 1:length(gs)){
    currgs=gs[[i]]

    currgsb = currgs$backward
    currgsf = currgs$forward
    names(currgsb)=get_DC_names(track,names(currgsb),i)
    names(currgsf)=get_DC_names(track,names(currgsf),i)
      
    back=lapply(names(currgsb),function(j){
      currgroup = currgsb[[j]]
      if(is.null(currgroup)){
        return (NULL)
      }
      names(currgroup)=get_DC_names(track,names(currgroup),i-1)
      data.frame(timestep=i,group1=j,direction=I("backward"),timestep2=i-1,group2=names(currgroup),similarity=unlist(currgroup))
    })
    back=do.call(rbind,back)
    
    forward=lapply(names(currgsf),function(j){
      currgroup = currgsf[[j]]
      if(is.null(currgroup)){
        return (NULL)
      }
      names(currgroup)=get_DC_names(track,names(currgroup),i+1)
      data.frame(timestep=i,group1=j,direction=I("forward"),timestep2=i+1,group2=names(currgroup),similarity=unlist(currgroup))
    })
    forward=do.call(rbind,forward)
    allsim=rbind(allsim,back,forward)
  }
  return (allsim)
}

get_flux_colors=function(track,allcols,cols2,dcmembership=NULL,singlecol=F,movecol="red",bysource=T,singlecolremain=T,remaincol="grey"){
  #get a per slice colour vector set up for python
  

  allflux=move_events_df(track,dcmembership,T)
  
  fluxcols1=lapply(unique(allflux$slice),function(x){
    currslice=allflux[allflux$slice==x,]
    
    currflux=data.frame(time=x-2,source=currslice$parent,target=currslice$child,fromlab=get_cluster_names(track,currslice$parent,x-1),tolab=get_cluster_names(track,currslice$child,x))
    if(singlecol){
      currflux$col=rgb(t(col2rgb(movecol)),maxColorValue = 255)
    }
    if(bysource&!singlecol){
      currflux$col=allcols[match(as.numeric(currslice$parent),track$comm_all)]
    }else if(!singlecol){
      currflux$col=allcols[match(as.numeric(currslice$child),track$comm_all)]      
    }
    if(singlecolremain){
      
    #remains in grey
      currflux$col[currslice$parent==currslice$child]=rgb(t(col2rgb(remaincol)),maxColorValue = 255)
    }
    currflux
  })
  fluxcols1=do.call(rbind,fluxcols1)
  return (fluxcols1)
}

coldictionary=function(track,allcols){
  cols1=lapply(1:length(track$dcs),function(x){
    currdcs=track$dcs[[x]]#
    cols2=allcols[match(currdcs,track$comm_all)]
    py_dict(cols2,keys=0:(length(currdcs)-1))
  })
  
  #nest this dictionary
  cols2=py_dict(cols1,keys=0:(length(track$dcs)-1))
  return(cols2)
}

get_alluvialplot=function(track,dcmembership,allcols,fluxbysource=T,fluxsinglecol=T,fluxmovecol="grey",fluxsinglecolremain=T,fluxremaincol="grey",
         fluxalpha=0.4,figwidth=8,figheight=2,rlabels=NULL,rstart=NULL,rstop=NULL,
         rmargins=c(0,0.2,1,1),
         cwidth=0.2,clusterlw=0.5,
         labelsize=0,
         reimport=T,removefile=T,exportfilename="Rplot.png")
  {
  
  cols2=coldictionary(track,allcols)

  fluxcols1=get_flux_colors(track=track,dcmembership=dcmembership,allcols,cols2,
                            fluxbysource,
                          singlecol=fluxsinglecol,movecol=fluxmovecol,
                          singlecolremain=fluxsinglecolremain,remaincol=fluxremaincol)

    if(is.null(rstop)){
      rstop=length(track$dcs)
    }
    if(is.null(rstart)){
      rstart=0
    }
    if(is.null(rlabels)){
      rlabels=c(1:length(track$dcs))
    }
	
    source_python("MTprocess.py")
    R_make_figure(track,cols2,figwidth,figheight,rmargins=rmargins,rstart=rstart,rstop=rstop,cwidth=cwidth,clusterlw=clusterlw,rlabels=rlabels,
                exportfilename=exportfilename,labelsize=labelsize,
                fluxalpha=fluxalpha,fluxfacecolor=fluxcols1$col,fluxfacefrom=fluxcols1$fromlab,fluxfaceto=fluxcols1$tolab,fluxfacets=fluxcols1$time)
    if(reimport){
      
      alluplot=load.image(exportfilename)
      par(mai=c(0,0,0,0))
      plot(alluplot,axes=F,rescale=T,xaxs="i",yaxs="i",asp="varying")
      rimsize=dev.size("in")
      
      par(mai=c(rmargins[2]*rimsize[2],
                (rmargins[1]*rimsize[1])+(cwidth*((1+rstart)-1)),
                (1-rmargins[4])*rimsize[2],
                ((1-rmargins[3])*rimsize[1])+(cwidth*((rstop)-length(track$dcs)))
                ),new=T)
      plot(0,type="n",xaxs="i",yaxs="i",xlim=c(rstart+1.5,rstop+1.5)-1,ylim=c(0,1),xlab="",ylab="",axes=F)
    }
    if(removefile&file.exists(exportfilename)){
      invisible(file.remove(exportfilename))
    }
  
  }


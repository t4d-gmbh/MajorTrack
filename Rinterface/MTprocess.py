from majortrack import MajorTrack
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from copy import deepcopy

def R_do_track(groupings,individuals,history):
	mt = MajorTrack(
		clusterings=groupings,
		individuals=individuals,
		history=history,
	)
	mt.get_group_matchup('fraction')
	# create the different instances with different history parameters

	mt.get_dcs()
	mt.get_community_membership()
	mt.get_community_group_membership()
	mt.get_individual_membership()
	mt.get_individual_group_membership()
	mt.get_community_events()
	mt.get_community_lifespans()
	mt.get_community_merges()
	mt.get_community_splits()

	return(mt)
	
def R_make_figure(mt,cols,figwidth,figheight,rmargins,rstop,rlabels,exportfilename,labelsize=0,rstart=0,l_size = 12,cwidth=0.2,clusteralpha=1,clusterlw=0.5,fluxalpha=0.4,fluxlw=0,
	fluxfacecolor='cluster',fluxfacefrom=None,fluxfaceto=None,fluxfacets=None):

	SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 12, 14, 18
	# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
	
	
	if type (fluxfacecolor) is list and fluxfacefrom is not None:
		#build dict
		new_cols=dict()
		for i in set(fluxfacets):
			i=int(i)
			indexes=[j for j in range(0,len(fluxfacets)) if fluxfacets[j]==i]
			temp_dict=dict()
			for j in indexes:
				fluxfacefrom[j]
				temp_dict[int(fluxfacefrom[j]),int(fluxfaceto[j])]=fluxfacecolor[j]
			new_cols[i]=temp_dict
		#print (new_cols)
		fluxfacecolor=new_cols
	

	#define figure size 
	fig1 = plt.figure(figsize=(figwidth,figheight),dpi=600)

	ax = fig1.add_subplot(1,1,1)
	ax.axis('equal')
	with_xaxis=True
	plt.margins(0)

	mt.get_alluvialdiagram(
		ax,
		invisible_x=not with_xaxis,
		cluster_width= cwidth,
		cluster_facecolor=cols,
		cluster_edgecolor=[0,0,0],
		with_cluster_labels= True,
		cluster_label= 'group_index',
		cluster_label_margin= (0, 1),
		x_axis_offset= 1,
		redistribute_vertically= 1,
		cluster_kwargs= {'alpha': clusteralpha, 'lw': clusterlw},
		label_kwargs= {'fontweight': 'light','fontsize': labelsize},
		flux_facecolor=fluxfacecolor,
		flux_kwargs= {'alpha': fluxalpha, 'lw': fluxlw}
		)

	ax.set_aspect('auto')
	ax.set_xlim([rstart, rstop])
	plt.subplots_adjust(left=rmargins[0], bottom=rmargins[1], right=rmargins[2], top=rmargins[3], wspace=0, hspace=0)
	rlabels=list(rlabels)

	if with_xaxis:
		tp = [(t + .5*(mt.slice_widths[i])) + 0*cwidth
			for i, t in enumerate(mt.timepoints)]
		ax.set_xticks(tp, minor=False)
		ax.set_xticklabels(rlabels,
			minor=False,
			size=l_size,
			horizontalalignment='center'
			)
		ax.tick_params(axis=u'x', which=u'both', length=0)
		plt.setp(ax.get_xticklabels(), visible=True)
	
	fig1.savefig(exportfilename)
	plt.close('all')

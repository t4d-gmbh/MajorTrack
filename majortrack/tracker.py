#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
"""
import itertools
from functools import wraps
from collections import Counter
from datetime import timedelta

from .lazylist import lazy_list
from .similarities import group_similarity_fraction, group_similarity_jaccard

from colorseq import DistinctColors
from pyalluv import AlluvialPlot, Cluster, Flux


class MajorTrack(object):
    r"""

    Parameters
    ===========

    clusterings: list, dict
      Sequence of clusterings. Each clustering must be present in the form of
      a list of membership-sets, i.e. a list of clusters with each cluster
      being defined by a `set` of data source associated to it (its members).

      **If provided as a `dict`**:

      keys: float, datetime
        The time points.
      values: list, dict
        The membership list of each clustering indicating to which cluster
        a data source belongs.
        See :obj:`~MajorTrack.clusterings` for details.
    history: int
      sets the number of time points (or slices) the algorithm can maximally
      go back in time to check for majority matches.
    \**kwargs optional parameter:
      timepoints: list
        The time points of each clustering.

        .. note::
          If `clusterings` if of type `dict` then the keys will be used as time
          points and this optional parameter is ignored, even if provided.

      slice_widths: list, float (default=None)
        The temporal duration of each snapshot in the sequence of clusterings.
        If not provided then simply the difference between time point `i` and
        `i+1` is used as the width of slice `i`. The width of the last slice is
        assumed to be the same as the duration of the 2nd last slice.
      individuals: list
        A list of all distinct data sources present in the dataset.

        .. todo::

          Build it from `self.clusterings`.

      group_matchup_method: str (default='fraction')
        Set the method to calculate the similarity between two clusters from
        different clusterings. By default the fraction of identical members is
        used as explained in the original article :cite:`liechti2019time`.
      use_lazylists: bool (default=False)
        Determine if :obj:`~LazyList`'s should be used to store data about
        dynamic clusters or normal lists.
        Most likely you want to use normal lists.

    .. bibliography:: ../references.bib

    Attributes
    ==========
    clusterings: list(list(set)) or LazyList
      Holds for each time point the configuration of the respective clustering.
      The clustering is given by a list of member-`set`s, with each `set`
      containing the data sources in a cluster.
    dcs: list or LazyList
      Ensemble of all dynamic clusters.

      .. todo::

          What's the type of an element? Is it just an identifier?
    length: int
      Number of slices present in the dataset.
    cluster_trace: list or LazyList
      Ensemble of all tracing paths the dynamic clusters.

      .. todo::

          What's the type of an element?
    group_matchup: list
      Holds for each time point the tracing and mapping sets of all clusters.
      Each element is a `dict` with the keys ``'forward'`` and ``'backward'``.
      Both hold a `dict` indicating for a cluster the best matching cluster
      along with the similarity score of the particular relation in a `tuple`.

      Example
      -------
      .. code-block:: python

         self.group_matchup[1] = {
             'backward': {0: (0, 1.0), ...},
         #                ^   ^  ^
         #                |   |  similarity score
         #                |   cluster from previous time point
         #                cluster from current time point.
             }

    group_mappings: list(list)
      Holds for each slice a list of mapping sets. The list is ordered like
      :obj:`~.MajorTrack.clusterings`.

      Example
      --------
      .. code-block:: python

          mt = MajorTrack(...)
          idx, cluster_id = 0, 1
          # get the set of data sources in this cluster
          c_members = mt.clusterings[0][1]
          # get the corresponding mapping set (from index idx + 1)
          mapping_set = mt.group_mappings[0][1]

    group_tracings: list(list)
      Holds for each slice a list of tracing sets. The list is ordered like
      :attr:`~.MajorTrack.clusterings`.
    group_mappers: list(list)
      Holds for each slice a list of mapper sets. The list is ordered like
      :attr:`~.MajorTrack.clusterings`.
    group_tracers: list(list)
      Holds for each slice a list of tracer sets. The list is ordered like
      :attr:`~.MajorTrack.clustering`.
    comm_nbr: int
      Number of dynamic clusters present in dataset.
    comm_all: list(dc index what type?)
      List of all dynamic clusters  present in the dataset.
    comm_group_members: ?

      .. todo::

        Unsure about this.

    comm_members: list(dict)
      holding for each slice of the dataset a dictionary indicating for each
      cluster (`key`) a list of data sources (`values`).

      .. todo::

        Rename to `dc_members`

    individuals: list
      holds all data sources.
    individual_group_membership: list(dict)
      holding for each slice of the dataset a dictionary indicating for a data
      source the cluster it belongs to.
    individual_membership: list(dict)
      holding for each slice of the dataset a dictionary indicating for a data
      source the dynamic cluster it belongs to.
    community_births: list(tuple)
      holding all dynamic cluster birth events.

      .. todo::

        Check and report format of this attribute.

    community_deaths: list(tuple)
      A list holding all dynamic cluster death events.

      .. todo::

        Check and report format of this attribute.

    community_lifespans: dict
      providing for each dynamic cluster the lifespan in the unit slices:
      ``{comm_id: nbr_slices}``
    community_splits: list(list)
      holds all split events of dynamic clusters.
    community_cby_splits: list(list)
      dynamic clusters that occurred through a split.
    community_cby_split_merges: list(list)
      dynamic clusters that occurred through a split-merge event.
    community_dby_splits: list(list)
      dynamic clusters that vanished after a split event.
    community_dby_split_merges: list(list)
      dynamic clusters that vanished after a split-merge event.
    community_merges: list(list)
      holds all merge events of dynamic clusters.
    community_cby_merges: list(list)
      dynamic clusters that occurred through a merge event.
    community_dby_merges: list(list)
      dynamic clusters that vanished after a merge event.
    community_growths: list(list)
      reports all growth events, i.e. changes in the size of a dynamic cluster
      that are not related to split or merge events.
    community_shrinkages: list(list)
      reports all shrinkage events, i.e. decreases in the size of a dynamic
      cluster that are not related to split or merge events.
    community_autocorrs: dict(dc index, list)
      hold for each dynamic cluster a dictionary with the auto-correlation
      (`value`) between the index of a slice (`key`) and the previous slice.
      The autocorrelation is given by:

      .. math::

        \frac{|dc_{i} \cap dc_{j}|}{|dc_{i} \cup dc_{j}|_{res}}

      where :math:`i, j` are the indices `from_dix` and `to_idx` and
      :math:`|<selection>|_{res}` is the number of data sources within
      `<selection>` counting all data sources (if `residents=False`) or
      only those present in both slices (`residents=True`).

      .. todo:

        Add ref to method (for details see...)

    """
    def __init__(self, clusterings, history, **kwargs):
        self.dcs = None
        self.comm_nbr = None
        # related to colouring
        self.color_sequence = None
        self.comm_colours = None
        self.sp_community_colour_idx = None
        self.history = history
        self._use_lazy = kwargs.get('use_lazylists', False)
        assert isinstance(clusterings, (list, dict))
        if isinstance(clusterings, list):
            self.timepoints = kwargs.pop(
                    'timepoints', list(range(len(clusterings)))
                    )
            # this is former groupings
            self.clusterings = clusterings
            # sort both clusterings and timepoints according to timepoints
            self.timepoints, self.clusterings = zip(
                    *sorted(
                        zip(self.timepoints, self.clusterings),
                        key=lambda x: x[0]
                        )
                    )
        else:
            self.timepoints = sorted(clusterings.keys())
            self.clusterings = list(clusterings[tp] for tp in self.timepoints)
        self.length = len(self.timepoints)
        # now determine the slice widths
        self.slice_widths = kwargs.get('slice_widths', None)
        if isinstance(self.slice_widths, (float, int, timedelta)):
            self.slice_widths = [self.slice_widths for idx in self]
        elif isinstance(self.slice_widths, (list, tuple)):
            self.slice_widths = list(self.slice_widths)
        else:
            self.slice_widths = [
                    self.timepoints[i] - self.timepoints[p]
                    for p, i in self._pair_iter()
                    ]
            self.slice_widths.append(self.slice_widths[-1])

        # todo: set individuals via set from self.clusterings of not provided.
        self.individuals = kwargs.get('individuals', None)

        self.group_matchup_method = kwargs.get(
                    'group_matchup_method',
                    'fraction'
                    )
        self._community_iter = None
        self.individual_membership = None

    def __iter__(self):
        return iter(range(self.length))

    def _pair_iter(self, func=None, *args, **kwargs):
        if func is None:
            def func(self, prev, item, *args, **kwargs): return prev, item
        _iter = iter(range(self.length))
        prev = next(_iter)
        for item in _iter:
            yield func(self, prev, item, *args, **kwargs)
            prev = item

    def _as_pair_iterator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            nbr_args = len(args)
            fargs = ()
            if not nbr_args:
                if not kwargs:
                    return self._pair_iter(func)
                else:
                    idx_prev = kwargs.pop('idx_prev', None)
                    idx_next = kwargs.pop('idx_next', None)
            else:
                idx_prev = args[0]
                if nbr_args > 1:
                    idx_next = args[1]
                    if nbr_args > 2:
                        fargs = args[2:]
                else:
                    idx_next = None
            if idx_prev is None:
                if idx_next is None:
                    return self._pair_iter(func, *fargs, **kwargs)
                    idx_next = None
                else:
                    idx_prev = idx_next - 1
                    assert idx_prev >= 0
            elif idx_next is None:
                idx_next = idx_prev + 1
                assert idx_next < self.length
            return func(self, idx_prev, idx_next, *fargs, **kwargs)
        return wrapper

    @property
    def _resident_population(self,):
        return sorted(list(
            self._prev_individuals.intersection(
                self._next_individuals
                )
            ))

    @property
    def _combined_population(self,):
        return sorted(list(
            self._prev_individuals.union(
                self._next_individuals
                )
            ))

    @_as_pair_iterator
    def resident_population(
            self, idx_prev=None, idx_next=None,
            *args, **kwargs):
        """
        Return the resident population between two time points.

        The resident population is simply the intersect of the populations at
        both time points.

        If not arguments are provided then an iterator is returned that gets
        the set of resident individuals between each slice.

        If only one index is provided then the other one will be completed,
        i.e. idx_prev = idx_next - 1 or idx_next = idx_prev + 1

        If further arguments are provided (all have to be unnamed), then the
        intersect is taken between all of these time points.

        Example::

          self.resident_population(2,4,5) will return the resident
          population between the time points 2, 4 and 5

        Parameters
        ==========
        idx_prev: int (default=None)
          index of one of the two data slices.
        idx_next: int (default=None)
          index of one of the two data slices.


        Returns
        =======
        resident_population: set
          contains all data sources that are in both data slices.
        """
        self._prev_individuals = self.individuals[idx_prev]
        self._next_individuals = self.individuals[idx_next]
        resident_population = set(self._resident_population)
        self._prev_individuals, self._next_individuals = None, None
        if args:
            for arg in args:
                resident_population = resident_population.intersection(
                        self.individuals[arg]
                        )
        return resident_population

    @_as_pair_iterator
    def combined_population(
            self, idx_prev=None, idx_next=None,
            *args, **kwargs):
        """
        Returns combination of the populations of two (or more) time points.

        This is simply the union of the populations at both time points.
        If not arguments are provided then an iterator is returned that gets
        the set of combined individuals between each slice.

        If only one index is provided then the other one will be completed,
        i.e. idx_prev = idx_next - 1 or idx_next = idx_prev + 1

        If further arguments are provided (all have to be unnamed), then the
        union is taken between all of these time points.

        Example
        -------
        .. code-block:: python

          self.resident_population(2,4,5)

        This will return the combined population of the time points 2, 4
        and 5.

        Parameters
        ==========
        idx_prev: int (default=None)
          index of the 1st time point to get the population from.
        idx_next: int (default=None)
          index of the 2nd time point to get the population from.

          .. note::

            If both `idx_prev` and `idx_next` are `None` then a pairwise
            iterator is returned that allows to loop over the combined
            population of neighbouring time points.

        """
        self._prev_individuals = self.individuals[idx_prev]
        self._next_individuals = self.individuals[idx_next]
        combined_population = set(self._combined_population)
        self._prev_individuals, self._next_individuals = None, None
        if args:
            for arg in args:
                combined_population = combined_population.union(
                        self.individuals[arg]
                        )
        return combined_population

    @_as_pair_iterator
    def resident_fraction(self, idx_prev=None, idx_next=None, *args):
        """
        Get the fraction of the combined population of tow (or more) slices.

        This indicates the population fraction that is present at all time
        points (or slices).

        This is simply the size of the intersection divided by the size of
        the union of the populations
        If further arguments are provided (all have to be unnamed), then the
        resident fraction is computed between all of these time points.

        Parameters
        ==========
        idx_pres: int (default=None)
          index of a slice.

          .. note::

            If no index is provided then an iterator is returned that yields
            the resident fraction of the data sources present in neighbouring
            time points.

        idx_next: int (default=None)
          index of a slice.

        Returns
        =======
        resident_fraction: float, iterator
          indicating the fraction of the population of data sources (union of
          all) that is present in all slices.
          If no values for the parameters `idx_prev` and `idx_next` are
          provided this method returns an iterator that will yield the fraction
          of the resident population between any two consecutive slices.
        """
        self._prev_individuals = self.individuals[idx_prev]
        self._next_individuals = self.individuals[idx_next]
        resident_population = set(self._resident_population)
        combined_population = set(self._combined_population)
        self._prev_individuals, self._next_individuals = None, None
        if args:
            for arg in args:
                resident_population = resident_population.intersection(
                        self.individuals[arg]
                        )
                combined_population = combined_population.union(
                        self.individuals[arg]
                        )
        return len(resident_population) / float(len(combined_population))

    def _set_dcs(self,):
        """
        Reset the list (or :ref:`~LazyList`) of dynamic clusters (`self.dcs`),
        as well as the list (or :ref:`~LazyList`) of all tracing paths for a
        cluster (self.cluster_trace)

        .. todo::

          Is this description of cluster_trace accurate?
        """
        if self._use_lazy:
            self.dcs = lazy_list(
                    lazy_name='dcs'
                    )
            self.cluster_trace = lazy_list(
                    lazy_name='cluster_trace'
                    )
        else:
            self.dcs = []
            self.cluster_trace = []
        self._community_iter = itertools.count(
                start=0,
                step=1
                )

    def _reset_container(self, container):
        """
        """
        s_container = self.__dict__.get(container, None)
        if s_container is not None:
            if isinstance(s_container, lazy_list):
                s_container = lazy_list(
                        lazy_name=container
                        )
            else:
                s_container = []

    def _grouping_similarity(
            self,
            method,
            mode='forward',
            use_union=False,
            **kwargs):
        """
        Get the similarity between clusters from current clusterings.

        The current clusterings are determined by `self._prev_grouping` and
        `self._next_grouping` using the specified method.

        Parameters
        ==========
        method: str
          determines the method upon which the similarity quantity is based.
        mode: str (default='forward')
          indicate the temporal direction of for the asymmetric similarity
          measure.
          .. todo::

            Specify what `'forward'` means for the denominator.

        use_union: bool (default=True)
          If `use_union==True` then if mode is 'forward' for each group of the
          previous grouping the groups from the next grouping are compared to
          the union of the previous group with all groups from the next
          grouping including members from the previous group.
          If mode is 'backward' then the same applies with inverted roles for
          prev and next.

        Returns
        =======
        similarities: dict(dict)
          for all clusters from the starting clustering (`key`) a dictionary
          (`value`) providing for all clusters from the ending clustering
          (inner `key`) as similarity value (inner `value`)
        """
        # filter for resident population only
        _residents = self._resident_population
        _prev_grouping_res = [
                {
                    _id for _id
                    in _group
                    if _id in _residents
                    }
                for _group in self._prev_grouping
                ]
        _next_grouping_res = [
                {
                    _id
                    for _id in _group
                    if _id in _residents
                    }
                for _group in self._next_grouping
        ]
        # set the correct order
        if mode == 'backward':
            from_grouping_all = self._next_grouping
            from_grouping = _next_grouping_res
            from_grouping_idx = range(len(self._next_grouping))
            to_grouping_all = self._prev_grouping
            to_grouping = _prev_grouping_res
            to_grouping_idx = range(len(self._prev_grouping))
        else:
            from_grouping_all = self._prev_grouping
            from_grouping = _prev_grouping_res
            from_grouping_idx = range(len(self._prev_grouping))
            to_grouping_all = self._next_grouping
            to_grouping = _next_grouping_res
            to_grouping_idx = range(len(self._next_grouping))
        # now compare the clusterings
        # initialize the similarities
        similarities = {
                i: {
                    j: None
                    for j in range(len(to_grouping_all))
                    }
                for i in range(len(from_grouping_all))
                }
        if method.lower() == 'fraction':
            for i in from_grouping_idx:
                similarities[i] = {}
                from_group = from_grouping[i]
                # calculate similarity between both reference groups/unions
                for j in to_grouping_idx:
                    to_group = to_grouping[j]
                    similarities[i][j] = group_similarity_fraction(
                            to_group,
                            from_group,
                            )
        elif method.lower() == 'jaccard':
            # TODO: Something is wrong here!
            for i in from_grouping_idx:
                similarities[i] = {}
                from_group = from_grouping[i]
                # no unions just 1-by-1
                if not use_union:
                    from_reference = from_group
                # union of the to groups
                else:
                    _from_union = set(from_group)
                    for j in to_grouping_idx:
                        # this is the index of the group:
                        _to_group = to_grouping[j]
                        if _from_union.intersection(_to_group):
                            _from_union = _from_union.union(_to_group)
                    from_reference = _from_union
                # calculate similarity between both reference groups/unions
                for j in to_grouping_idx:
                    to_group = to_grouping[j]
                    similarities[i][j] = group_similarity_jaccard(
                            from_group,
                            to_group,
                            len(from_reference)
                            )
        else:
            # TODO: Implement mutual information (and others?)
            pass
        return similarities

    @_as_pair_iterator
    def _get_group_similarities(
            self, idx_prev=None, idx_next=None,
            **kwargs
            ):
        r"""
        A return similarities between the clusterings of two time points.

        Approach:
            1. Forward mapping:
                For each group from the previous time point, get the most
                similar group from the next time point. If there is no group in
                the next time point similar to a previous group, the group gets
                None associated.
            2. Backward tracing:
                For each group of the next time point get the most similar
                group from the previous time point.
            3. Match-up:
                - Both mapping and tracing agree: Match the groups
                - Tracing suggests a previous group that has a different group
                    as mapping: Probation mapping from next to previous group.
                    Probation mapping:
                        - If prev group gets no mapping to next
                            group, the probation mapping is considered valid.
                        - A probation mapping leads to a new group unless in
                            any next time step this group merges with the most
                            mapping (descendant) group of the previous group.
                        - Only the last probation mapping is kept and
                            considered when merging.

        Parameters
        ==========
        idx_prev: int (default=None)
          the previous (earlier) time stamp (slice).
        idx_next: int (default=None)
          the next (later) time stamp (slice)
        \**kwargs optional parameter:
          timepoints: list
          method: str (default='fraction')
            Set the method to calculate the similarity between two groups from
            different clusterings.
            Default option is the fraction of data sources
          min_overlap: int (default=1)
            sets the minimal number of data sources that two clusters need to
            share for them to be matched up (if the similarity is maximal).
          single_next: bool (default=False)
            If set to True then similarity is compared 1-by-1
            between a group from the previous time stamp and any groups
            from the next time stamp for step 1.
            If set to False (default) then the similarity is compared
            between groups from the next time stamp and the union of the
            previous group and all groups from the next time stamp
            containing a member from the previous group.
          single_prev: bool (default=False)
            If set to True then similarity is compared 1-by-1
            between a group from the next time stamp and any groups
            from the previous time stamp. If set to False (default) then
            the similarity is compared between groups from the previous
            time stamp and the union of the next group and all groups from
            the previous time stamp containing a member from the previous
            group.

        Returns
        =======
        similarities: dict(str, dict)
          holding two keys (`'forward'` and `'backward'`) with each holding a
          dictionary with the best matches between clusters
          (inner `key`: `'matchup'`) and the similarities
          (inner `key`: `'similarities'` )

        """
        method = kwargs.get('method', 'fraction')
        min_overlap = kwargs.get('min_overlap', 1)
        _union_next = not kwargs.get('single_next', False)
        _union_prev = not kwargs.get('single_prev', False)
        self._prev_grouping = self.clusterings[idx_prev]
        self._next_grouping = self.clusterings[idx_next]
        # todo: this is still from the TimeSeriesData
        # else:
        #     self._prev_graph = self.graphs[idx_prev]
        #     self._next_graph = self.graphs[idx_next]
        #     self._prev_grouping = self._get_grouping(self._prev_graph)
        #     self._next_grouping = self._get_grouping(self._next_graph)
        self._prev_individuals = self.individuals[idx_prev]
        self._next_individuals = self.individuals[idx_next]
        # Step 1: Forward mapping
        fw_similarities = self._grouping_similarity(
                method=method,
                mode='forward',
                use_union=_union_next
                )
        # Step 2: Backward tracing
        bw_similarities = self._grouping_similarity(
                method=method,
                mode='backward',
                use_union=_union_prev
                )
        # - Get most similar group mapping
        fw_matchup = self._matchup(
            fw_similarities,
            self._prev_grouping,
            self._next_grouping,
            min_overlap=min_overlap
        )
        # - Get most similar group tracing
        bw_matchup = self._matchup(
            bw_similarities,
            self._next_grouping,
            self._prev_grouping,
            min_overlap=min_overlap
        )
        self._prev_graph, self._next_graph = None, None
        self._prev_grouping, self._next_grouping = None, None
        self._prev_individuals, self._next_individuals = None, None
        return {
                'forward': {
                    'matchup': fw_matchup,
                    'similarities': fw_similarities
                    },
                'backward': {
                    'matchup': bw_matchup,
                    'similarities': bw_similarities
                    }
                }

    @staticmethod
    def _matchup(similarities, from_grouping, to_grouping, min_overlap):
        """
        Find the best match for all clusters between snapshots.

        This matches all clusters from `from_grouping` with clusters form
        `to_grouping`.
        If no match is found or the overlap in number of
        data sources is smaller than `min_overlap`, the match is set to None

        Parameters
        ==========
        similarities:
          .. todo::

            Description

        from_grouping:
          .. todo::

            Description

        to_grouping:
          .. todo::

            Description

        min_overlap:
          Set `min_overlap=0` to ignore the minimal overlap size.
          .. todo::

            Description


        Returns
        =======
        matchup: dict
          .. todo::

            Description

        """
        matchup = {}
        for i in range(len(from_grouping)):
            max_sim = max(
                    similarities[i].items(),
                    key=lambda x: x[1] if x[1] is not None else 0
                    )
            matchup[i] = ([max_sim[0]], max_sim[1])
            # assert the minimal overlap
            if min_overlap:
                ok_matchups = []
                for mg in matchup[i][0]:
                    if len(
                            from_grouping[i].intersection(
                                to_grouping[mg]
                                )
                            ) >= min_overlap:
                        ok_matchups.append(mg)
                matchup[i] = (
                        ok_matchups,
                        matchup[i][1] if ok_matchups else None
                        )
            # without min_overlap, just assert that max similarity > 0
            else:
                if matchup[i][1] == 0.0:
                    matchup[i] = ([], None)
            # check how many groups from the next timestep have the same value
            _sim_counts = {
                    sim: count
                    for sim, count in Counter(
                        similarities[i].values()
                        ).items()
                    }
            # there is more than 1 from next with a maximal similarity
            if matchup[i][0] and _sim_counts[matchup[i][1]] > 1:
                _other_next_groups = [
                    _gid
                    for _gid, _sim
                    in similarities[i].items()
                    if _sim == matchup[i][1]
                    and _gid not in matchup[i][0]
                ]
                # make the match-up a list of ids if there are several.
                if _other_next_groups:
                    matchup[i][0].extend(_other_next_groups)
        return matchup

    def get_group_matchup(self, matchup_method=None):
        r"""
        Determine majority relation between neighbouring snapshots.

        Parameters
        ===========
        matchup_method: str (default=None)
          If provided this overwrites
          :attr:`~majortrack.MajorTrack.group_matchup_method`.
          It determines the method to use when calculating similarities between
          clusters from neighbouring snapshots.

        Returns
        =======
        self: :class:`.MajorTrack`
          with new attribute :obj:`~.MajorTrack.group_matchup`.

        """
        if matchup_method is None:
            matchup_method = self.group_matchup_method
        # if self.group_matchup:
        self.group_matchup = []
        # if self.group_similarities:
        self.group_similarities = []
        self.group_matchup.append(
                {
                    'backward': {
                        _group_id: (None, None)
                        for _group_id in range(len(self.clusterings[0]))
                        }
                    }
                )
        self.group_similarities.append(
                {
                    'backward': {
                        _group_id: None
                        for _group_id in range(len(self.clusterings[0]))
                        }
                    }
                )
        for _idx in range(self.length - 1):
            _group_similarities = self._get_group_similarities(
                    _idx, _idx + 1,
                    method=matchup_method
                    )
            # set forward matchup/similarities for current step
            self.group_matchup[-1][
                    'forward'
                    ] = _group_similarities['forward']['matchup']
            self.group_similarities[-1][
                    'forward'
                    ] = _group_similarities['forward']['similarities']
            # create backward matchup/similarities for next step
            self.group_matchup.append(
                    {'backward': _group_similarities[
                        'backward'
                        ]['matchup']}
                    )
            self.group_similarities.append(
                    {'backward': _group_similarities[
                        'backward'
                        ]['similarities']}
                    )
        # complete forward matchup/similarites with None's
        self.group_matchup[-1]['forward'] = {
                _group_id: (None, None)
                for _group_id in range(len(self.clusterings[-1]))
                }
        self.group_similarities[-1]['forward'] = {
                _group_id: None
                for _group_id in range(len(self.clusterings[-1]))
                }

    # todo: method unused...
    def get_span(self, idx, span_set, get_indivs=True):
        r"""
        Create the tracer tree.

        Parameters
        ===========
        idx: int
          index of the slice in which to start.
        span_set: int, str
          If an `int` is provided it specifies the index of the target cluster.
          If a `str` is given, it is considered as the label of a data source
          and the containing cluster is selected.

          .. todo::

            The label of a cluster should be the only option.

        get_indivs: bool (default=True)
          If set to `True` a list of sets of individual is returned for each
          slice starting from the index.
          If it is set to `False` a list of cluster labels is returned for each
          slice.

        """
        span_tree = {}
        if isinstance(span_set, int):
            span_tree[idx] = [self.clusterings[idx][span_set]]
        elif isinstance(span_set, str):
            span_tree[idx] = filter(
                    lambda g: span_set in g,
                    self.clusterings[idx]
                    )
        else:
            span_tree[idx] = [span_set]
        current_set = set.union(*span_tree[idx])
        for _idx in range(idx + 1, self.length):
            next_groupings = self.clusterings[_idx]
            next_contained = [*filter(
                    lambda grp: any([memb in current_set for memb in grp]),
                    next_groupings
                    )]
            if next_contained:
                span_tree[_idx] = next_contained
                current_set = set.union(*next_contained)
            else:
                break
        if get_indivs:
            return span_tree
        else:
            span_tree_idxs = {}
            for _idx in span_tree:
                _span_set = set.union(*span_tree[_idx])
                span_tree_idxs[_idx] = [
                        i for i in range(len(self.clusterings[_idx]))
                        if any([
                            el in _span_set
                            for el in self.clusterings[_idx][i]
                            ])
                        ]
            return span_tree_idxs

    def get_dcs(self, bijective_paths=True, **kwargs):
        r"""
        Derives from the history of dynamic clusters from
        :attr:`~.MajorTrack.group_matchup`.

        .. todo::

          Rename to `get_dc`

        Parameters
        ===========
        bijective_paths: bool (default=True)
          If set to `True` then at each step in the construction of the tracing
          flow a mapping flow needs to map forward to the target cluster in
          order to continue to extend the tracing flow.

        \**kwargs optional parameter:
          from_idx: int
            starting index.

            .. note::

              At the starting index all clusters are per definition new
              dynamic clusters.

          to_idx: int
            Stopping index. The community detection algorithm will stop at
            this index (including it).
        """
        if self.dcs is None:
            self._set_dcs()
        else:
            self._reset_container('dcs')
            self._reset_container('cluster_trace')
        # Get the tracing/mapping tracer/mapper set
        self._get_tracings()
        self._get_mappings()
        self._get_mappers()
        self._get_tracers()
        if bijective_paths:
            def validate_path(*args):
                return self._from_flow(*args)
        else:
            def validate_path(*args):
                return True
        from_idx = kwargs.get('from_idx', 0)
        to_idx = kwargs.get('to_idx', self.length - 1)
        # TODO: from_idx must be 0 as otherwise the method might fail by
        # accessing an index < from_ids in self.dcs

        # populate communites with empty lists up to the starting point
        for _idx in range(from_idx):
            self.dcs.append([])
            self.cluster_trace.append([])
        for idx in range(from_idx, to_idx + 1):
            self.dcs.append([])
            self.cluster_trace.append([])
            for group_id in range(len(self.clusterings[idx])):
                s_c = None
                tracing_flow = self.get_flow(
                    idx, {group_id},
                    validate_path=validate_path
                )
                go_back = len(tracing_flow) - 1
                fwd_flow = []
                while go_back >= 1:
                    # get the potential source set
                    pot_source = tracing_flow[go_back]
                    # check if pot_source set has single community identity
                    if len(set(
                            (
                                self.dcs[idx-go_back][psg]
                                for psg in pot_source
                                )
                            )) == 1:
                        pot_mapping_flow = self.get_flow(
                                idx-go_back, pot_source, False, go_back,
                                validate_path=validate_path
                                )
                        if len(pot_mapping_flow) == go_back + 1:
                            if pot_mapping_flow[-1] == {group_id}:
                                s_c = self.dcs[idx-go_back][next(
                                    iter(pot_source)
                                    )]
                                fwd_flow = pot_mapping_flow
                                break
                    go_back -= 1
                if s_c is None:
                    s_c = next(self._community_iter)
                self.dcs[idx].append(s_c)
                self.cluster_trace[idx].append(0)
                # combine tracing_flow and forward flow > included_flow
                # relabel the included_flow
                included_flow = self._create_flow(
                        [tracing_flow[:len(fwd_flow)], fwd_flow[::-1]]
                        )
                if included_flow:
                    self._relabel_included_flow(
                            idx, s_c, included_flow
                            )
                # # now the marginal dcs still need to be relabeled
                marginal_flows = self.get_marginal_flows(idx, included_flow)
                self._relabel_marginal_flows(s_c, idx, marginal_flows)

    def _relabel_included_flow(self, idx, source_community, included_flow):
        r"""
        Associate a new dynamic cluster label to set of clusters in a slice.

        Parameters
        -----------
        idx: int
          index of the slice in which to relabel a dynamic cluster.
        source_community: str
          label of the dynamic cluster to use as new label.
        included_flow: list
          Ensemble of clusters from this slice to relabel.
        """
        incl_tree_length = len(included_flow)
        if incl_tree_length > 2:
            for step_back, inc_groups in enumerate(included_flow):
                if 0 < step_back < incl_tree_length - 1:
                    for inc_grp in inc_groups:
                        if self.dcs[
                                idx - step_back][
                                        inc_grp] != source_community:
                            self.dcs[
                                    idx - step_back
                                    ][inc_grp] = source_community
                            self.cluster_trace[idx - step_back][inc_grp] = 1

    # TODO: should we also use non-bijective paths?
    #       as `validate_path` is not provided here and thus defaults to
    #       self._from_flow (so bijective)
    def get_marginal_flows(self, idx, included_flows):
        r"""
        Determines the ensemble of marginal clusters given a target cluster and
        its identity flow.

        Parameters
        -----------
        idx: int
          index of the slice in which target cluster is situated.
        included_flows: list
          Identity flow of the target cluster.

        """
        flow_l = len(included_flows)
        if flow_l:
            source_set = included_flows[-1]
            target_set = included_flows[0]
            source_flow = self.get_flow(
                idx-flow_l+1, source_set, False, flow_l-1, majority=False,
                    )
            sl = len(source_flow)
            source_flow += [set() for _ in range(flow_l-sl)]
            target_flow = self.get_flow(
                idx, target_set, True, flow_l-1, majority=False,
                    )
            tl = len(target_flow)
            target_flow += [set() for _ in range(flow_l-tl)]
            return [
                    set.intersection(
                        source_flow[::-1][i], target_flow[i]
                        ).difference(included_flows[i])
                    for i in range(flow_l)
                ]
        else:
            return []
        # # Using step-wise tracers/mappers (old way)
        # marg_sources = list(included_flows)
        # marg_sinks = list(included_flows)
        # # get all potential sources, i.e. included_flows and all tracers
        # for _s in range(flow_l - 1):
        #     # set all marginal source
        #     c_source = marg_sources[flow_l - _s - 1]
        #     _idx = idx - flow_l + 1 + _s
        #     for csg in c_source:
        #         marg_sources[flow_l-_s-2] = marg_sources[flow_l-_s-2].union(
        #                 set(self.group_tracers[_idx][csg])
        #                 )
        #     c_sink = marg_sinks[_s]
        #     _idx = idx - _s
        #     for csg in c_sink:
        #         marg_sinks[_s+1] = marg_sinks[_s+1].union(
        #                 set(self.group_mappers[_idx][csg])
        #                 )
        # # only keep groups that are both in marg_sources and marg_sinks but
        # # not in the included flow
        # return [set.intersection(
        #     marg_sources[i], marg_sinks[i]).difference(included_flows[i])
        #     for i in range(flow_l)
        #     ]

    def _relabel_marginal_flows(self, source_community, idx, flows):
        r"""
        Associate all clusters in a marginal flow to a dynamic cluster.

        Parameters
        -----------
        source_community: str
          Label of the dynamic cluster   :ref:`flows`
          to the dynamic cluster
        idx: int
          index of the slice in which target cluster is situated.
        """
        for _i, marg_set in enumerate(flows):
            for mg in marg_set:
                self.dcs[idx-_i][mg] = source_community
                self.cluster_trace[idx-_i][mg] = -1

    def _from_flow(self, idx, source_set, flow, bwd=True):
        r"""
        Check if the source set stems form the provided tracing flow.

        Parameters
        -----------
        idx: int
          index of the slice in which the flow starts.
        source_set: set
          Set of clusters to test whether it is included in the flow or not.

          .. todo::

              Clarify the type of each element.
              Are they simply cluster labels or membership lists?

              How do we know at which index `source_set` is situated? Simply
              ``idx - len(flow)``?

        flow: list(set)
          Flow (either a tracing or a mapping flow) starting with the target
          cluster.
        bwd: bool (default=True)
          Indicating the temporal direction of the flow. If `bwd=True`
          the direction is `backwards` in time and the flow thus
          corresponds to the tracing flow. If set to `False` the flow
          is considered as a mapping flow.

        Returns
        -------
        is_from_flow: bool
          indicate whether a source set stemps from the provided flow or not.
        """
        reach = len(flow)
        if bwd:
            _get_next = self._get_tracing_sets
            ts_lim = 0
            _dir = -1
        else:
            _get_next = self._get_mapping_sets
            ts_lim = self.length
            _dir = 1
        step = 0
        _idx = idx + _dir * step
        next_set = source_set
        while step < reach and _dir * (ts_lim - _idx) > 0:
            next_set = _get_next(_idx, next_set)
            if next_set and next_set.issubset(flow[-1 - step]):
                return True
            step += 1
            _idx = idx + _dir * step
        return False

    def get_flow(self, idx, source_set, bwd=True, max_dist=None, **kwargs):
        r"""
        Parameters
        ----------
        idx: int
          time series index defining the starting point
        source_set: set
          set of clusters at the starting point slice.
        bwd: bool (default=True)
          indicating the direction, `True` is backward, `False` forward.
        max_dist: int (default=None)
          set the maximal length of the flow.
        \**kwargs optional parameter:
          majority: bool (default=True)
            allows to specify if of only the majority should be used to move
            between time-points.
          validate_path: function (default=:meth:`~.MajorTrack._from_flow`
            Provide a validation method to use during the construction of a
            flow.

        Returns
        -------
        flow: list
          identity flow starting (including) from the source set.
        """
        _maj = kwargs.get('majority', True)
        validate_path = kwargs.get('validate_path', self._from_flow)
        if bwd:
            if _maj:
                _get_next = self._get_tracing_sets
            else:
                _get_next = self._get_mappers_sets
            ts_lim = 0
            _dir = -1
        else:
            if _maj:
                _get_next = self._get_mapping_sets
            else:
                _get_next = self._get_tracers_sets
            ts_lim = self.length
            _dir = 1
        if max_dist is None:
            max_dist = self.history
        step = 0
        _idx = idx + _dir * step
        _new_set = source_set
        flow = [source_set]
        while step < max_dist and _dir * (ts_lim - _idx) > 0:
            # get for each group the next set
            _new_set = _get_next(_idx, _new_set)
            # only keep sets that stem from flow
            if _new_set and validate_path(_idx+_dir, _new_set, flow, not bwd):
                flow.append(_new_set)
                step += 1
                _idx = idx + _dir * step
            else:
                break
        return flow

    def _get_mapping_sets(self, idx, source_set):
        r"""
        Get the mapping set of the provided source set.

        .. todo::

            This is more an itemgetter of :attr:`~.MajorTrack.group_mappings`.

            Actually not really.

        Parameters
        -----------
        idx: int
          time series index defining the starting point
        source_set: set
          set of clusters at the starting point slice.

        Returns
        -------
        mapping_set: set
          the mapping set of the provided source set.
        """
        if source_set:
            return set.union(
                    *[set(self.group_mappings[idx][sg]) for sg in source_set]
                    )
        else:
            return set()

    def _get_tracing_sets(self, idx, source_set):
        r"""
        Get the tracing set of the provided source set.

        Parameters
        -----------
        idx: int
          time series index defining the starting point
        source_set: set
          set of clusters at the starting point slice.

        Returns
        -------
        tracing_set: set
          the tracing set of the provided source set.

        """
        return set.union(
                *[set(self.group_tracings[idx][sg]) for sg in source_set]
                )

    def _get_mappers_sets(self, idx, source_set):
        r"""
        idx: int
          time series index defining the starting point

        Parameters
        -----------
        idx: int
          time series index defining the starting point
        source_set: set
          set of clusters at the starting point slice.

        Returns
        -------
        mapper_set: set
          the mapper set of the provided source set.
        """
        return set.union(
                *[set(self.group_mappers[idx][sg]) for sg in source_set]
                )

    def _get_tracers_sets(self, idx, source_set):
        r"""

        Parameters
        -----------
        idx: int
          time series index defining the starting point
        source_set: set
          set of clusters at the starting point slice.

        Returns
        -------
        tracer_set: set
          the tracer set of the provided source set.
        """
        return set.union(
                *[set(self.group_tracers[idx][sg]) for sg in source_set]
                )

    def _get_tracings(self,):
        r"""
        Compute the tracing sets for all clusters in the sequence of slices.

        Returns
        --------
        None: NoneType
          Set new attribute :ref:`group_tracings`.

        """
        self.group_tracings = []
        for idx in range(self.length):
            if idx:
                _current_bw_matchup = self.group_matchup[idx]['backward']
                self.group_tracings.append(
                    [
                        _current_bw_matchup[_group][0]
                        for _group in range(len(self.clusterings[idx]))
                        ]
                )
            else:
                self.group_tracings.append(
                    [[] for _group in range(len(self.clusterings[idx]))]
                )

    def _get_mappings(self,):
        r"""
        Compute the mapping sets for all clusters in the sequence of slices.

        Returns
        --------
        None: NoneType
          Set new attribute :ref:`group_tracings`.

        """
        self.group_mappings = []
        for idx in range(self.length - 1):
            _current_fw_matchup = self.group_matchup[idx]['forward']
            self.group_mappings.append(
                [
                    _current_fw_matchup[_group][0]
                    for _group in range(len(self.clusterings[idx]))
                    ]
            )
        self.group_mappings.append(
            [[] for _group in range(len(self.clusterings[-1]))]
        )

    def _get_mappers(self,):
        r"""
        Compute the mapper sets for all clusters in the sequence of slices.

        Returns
        --------
        None: NoneType
          Set new attribute :ref:`group_mappers`.

        """
        self.group_mappers = [
            [[] for _group in range(len(self.clusterings[0]))]
        ]
        for idx in range(1, self.length):
            self.group_mappers.append(
                [[] for _group in range(len(self.clusterings[idx]))]
            )
            for prev_grp in range(len(self.group_mappings[idx - 1])):
                for mapps in self.group_mappings[idx - 1][prev_grp]:
                    self.group_mappers[idx][mapps].append(prev_grp)

    def _get_tracers(self,):
        r"""
        Compute the tracer sets for all clusters in the sequence of slices.

        Returns
        --------
        None: NoneType
          Set new attribute :ref:`group_tracers`.

        """
        self.group_tracers = []
        for idx in range(self.length - 1):
            self.group_tracers.append(
                [[] for _group in range(len(self.clusterings[idx]))]
            )
            for next_grp in range(len(self.group_tracings[idx + 1])):
                for tracing in self.group_tracings[idx + 1][next_grp]:
                    self.group_tracers[idx][tracing].append(next_grp)
        self.group_tracers.append([
            [[] for _group in range(len(self.clusterings[-1]))]
        ])

    def _create_flow(self, paths, contained=True):
        """
        Combine arbitrary many paths to a flow.

        Example
        --------

        .. code-block::

            [[{1}, {2}, {1}], [{0}, {1}, {1}]] > [{0, 1}, {1, 2}, {1}]

        Parameters
        -----------
        paths: list(list)
          List of individual paths to combine to a flow.
        contained: bool (default=True)
          If a single path contains a None at a specific point, then
          the flow is cut at this point.
          If contained is set do False the flow continues until the
          longest path.

        Returns
        --------
        flow: list
          The flow resulting from the list of paths provided.
        """
        if contained:
            flow = [
                    grps
                    for grps in map(
                            lambda x: set.union(*[*map(set, x)]),
                            zip(*paths)
                            )
                    if None not in grps
                    ]
        else:
            flow = [
                    grps
                    for grps in map(
                            lambda x: set.union(*[*map(set, x)]),
                            zip(*paths)
                            )
                    ]
            for el in flow:
                if None in el:
                    el.remove(None)
        return flow

    def get_community_group_membership(self,):
        """
        Defines per timepoint a list of clusters belonging to a dynamic cluster

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.comm_group_members`
          - :attr:`~.MajorTrack.comm_all`
          - :attr:`~.MajorTrack.comm_nbr`

        """
        _all_comms = []
        for idx in range(self.length):
            _all_comms.extend(
                    [com for com in self.dcs[idx]]
                        )
        self.comm_all = sorted(set(_all_comms))
        self.comm_nbr = len(self.comm_all)
        # get community members
        self.comm_group_members = []
        for idx in range(self.length):
            current_comms = {}
            for group_id in range(len(self.clusterings[idx])):
                try:
                    current_comms[self.dcs[idx][group_id]].append(
                            group_id
                            )
                except KeyError:
                    current_comms[self.dcs[idx][group_id]] = [group_id]
            self.comm_group_members.append(
                    current_comms
                    )

    def get_community_membership(self,):
        """
        Defines for each time point a membership list of data sources for each
        existing dynamic cluster.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.comm_members`
        """
        self.comm_members = []
        for idx in range(self.length):
            # todo: switch to defaultdict to avoid managing the exception
            current_comms = {}
            for group_id in range(len(self.clusterings[idx])):
                try:
                    current_comms[self.dcs[idx][group_id]].extend(
                            self.clusterings[idx][group_id]
                            )
                except KeyError:
                    current_comms[
                            self.dcs[idx][group_id]
                            ] = list(self.clusterings[idx][group_id])
            self.comm_members.append(
                    current_comms
                    )

    def get_individual_group_membership(self,):
        """
        Defines for each time point a dict holding for each data source its
        cluster membership.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.individual_group_membership`
        """
        self.individual_group_membership = []
        for idx in range(self.length):
            self.individual_group_membership.append(
                {
                    indiv: grp_id
                    for grp_id, members in enumerate(
                        self.clusterings[idx]
                        )
                    for indiv in members
                }
            )

    def get_individual_membership(self,):
        """
        Defines for each time point a dict holding for each data source its
        dynamic cluster membership.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.individual_membership`
        """
        self.individual_membership = []
        for idx in range(self.length):
            self.individual_membership.append(
                {
                    indiv: comm_id
                    for comm_id, members in self.comm_members[idx].items()
                    for indiv in members
                }
            )

    def get_community_births(self):
        """
        Determines all birth events.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_births`
        """
        self.community_births = []
        _births = []
        # consider all clusters at t0 as new births
        for com in set(self.dcs[0]):
            _births.append(([], [com]))
        self.community_births.append(_births)
        # check if there exist clusters with exclusively new members
        for i in range(1, self.length):
            _births = []
            res_pop = self.resident_population(i-1, i)
            for comm in self.comm_members[i]:
                if not any(
                    [memb in res_pop for memb in self.comm_members[i][comm]]
                ):
                    _births.append(
                        ([], [comm])
                    )
            self.community_births.append(_births)
        return None

    def get_community_deaths(self,):
        """
        Determines all dynamic community death events.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_deaths`

        """
        # at t0 there can be no deaths
        self.community_deaths = [[]]
        for i in range(1, self.length):
            _deaths = []
            res_pop = self.resident_population(i-1, i)
            for comm in self.comm_members[i-1]:
                # if there are no members in the next time step still present
                if not any(
                    [memb in res_pop for memb in self.comm_members[i-1][comm]]
                ):
                    _deaths.append(
                        ([comm], [])
                    )
            self.community_deaths.append(_deaths)
        return None

    def get_community_lifespans(self,):
        """
        Determines the lifespans of all dynamic clusters.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_lifespans`

        """
        _life_spans = {
            c: 0 for c in self.comm_all
        }
        for idx in range(self.length):
            for comm in set(self.dcs[idx]):
                _life_spans[comm] += 1
        self.community_lifespans = _life_spans
        return None

    def get_community_avg_lifespan(self, mode='ensemble'):
        """
        Determines the lifespans of all dynamic clusters.

        Parameters
        ===========
        mode: str (default='ensemble')
          Determines what type of average should be computed.
          Possible are either `ensemble` (default) or
          'weighted_per_indiv_per_slice'.
          The `ensemble` average simply consists of the arithmetic mean of all
          lifespans.
          The `weighted_per_indiv_per_slice` yields the average value of the
          life span of a dynamic cluster a randomly picked data source belongs
          to during at randomly picked slice.

        Returns
        =======
        avg_dc_lifespan: float
          the average number of slices a dynamic cluster exists.
        """
        self.get_community_lifespans()
        if mode == 'weighted_per_indiv_per_slice':
            # hold sum over size over slice
            _tot_counts = {c: 0 for c in self.comm_all}
            for idx in range(self.length):
                c_comms = set(self.dcs[idx])
                for c_com in c_comms:
                    # how many individuals are in currently
                    c_size = len(self.comm_members[idx][c_com])
                    _tot_counts[c_com] += c_size
            _grand_total = float(sum(_tot_counts.values()))
            _c_weight = {c: _tot_counts[c]/_grand_total for c in self.comm_all}
            return sum([
                _c_weight[c]*cls
                for c, cls in self.community_lifespans.items()
                ])
        else:
            return sum(
                    self.community_lifespans.values()
                    ) / float(self.comm_nbr)

    def get_community_splits(self,):
        """
        Get all split events and determine what clusters arise through a
        pure split event, i.e. not a split-merge combination.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_splits`
          - :attr:`~.MajorTrack.community_cby_splits`
          - :attr:`~.MajorTrack.community_cby_split_merges`
          - :attr:`~.MajorTrack.community_dby_splits`
          - :attr:`~.MajorTrack.community_dby_split_merges`

        """
        # no splits possible at the first time step
        self.community_splits = [[]]
        self.community_cby_splits = [[]]
        self.community_cby_split_merges = [[]]
        self.community_dby_splits = [[]]
        self.community_dby_split_merges = [[]]
        for i in range(1, self.length):
            _splits = []
            _new_cby_split = []
            _destroyed_by_split = []
            _new_cby_split_merge = []
            _destroyed_by_split_merge = []
            comms = [*self.comm_members[i-1].keys()]
            res_pop = self.resident_population(i-1, i)
            for comm, membs in self.comm_members[i-1].items():
                next_comms = set(
                        self.individual_membership[i].get(memb, None)
                        for memb in membs
                )
                if None in next_comms:
                    next_comms.remove(None)
                # if we find the members of `comm` at `i-1` in more than 1
                # community at `i`, we have a split.
                if len(next_comms) > 1:
                    _splits.append(
                        ([comm], list(next_comms))
                    )
                other_next_comms = [nc for nc in next_comms if nc != comm]
                if other_next_comms:
                    # check if any of the communities at `i` is new
                    new_next_comms = [
                            onc
                            for onc in other_next_comms
                            if onc not in comms
                            ]
                    # check if only members from the `comm` are in the new
                    # community/ies
                    for nnc in new_next_comms:
                        new_mc_res = [
                                nmc
                                for nmc in self.comm_members[i][nnc]
                                if nmc in res_pop
                                ]
                        if all([nmcr in membs for nmcr in new_mc_res]):
                            # new community through split from old community
                            _new_cby_split.append(
                                    ([comm], [nnc])
                                    )
                        else:
                            # the new community/ies at `i` also contain
                            # resident elements that did not belong to `comm`
                            # at `i-1`, thus communities from `i-1` have merged
                            # into the new community/ies.
                            orig_comms = set(
                                    self.individual_membership[i-1].get(nmcr)
                                    for nmcr in new_mc_res
                                    )
                            _new_cby_split_merge.append(
                                    (list(orig_comms), [nnc])
                                    )
                    # check if comm was destroyed
                    if comm not in next_comms:
                        # do all members of the new community come from `comm`?
                        next_membs = [
                                memb
                                for nc in next_comms
                                for memb in self.comm_members[i][nc]
                                if memb in res_pop
                                ]
                        if all([nm in membs for nm in next_membs]):
                            # all come from `comm` > split destruction
                            _destroyed_by_split.append(
                                    ([comm], [])
                                    )
                        else:
                            # also a merge > split-merge destruction
                            _destroyed_by_split_merge.append(
                                    ([comm], [])
                                    )
            self.community_cby_splits.append(_new_cby_split)
            self.community_dby_splits.append(_destroyed_by_split)
            self.community_cby_split_merges.append(_new_cby_split_merge)
            self.community_dby_split_merges.append(_destroyed_by_split_merge)
            self.community_splits.append(_splits)

    def get_community_merges(self,):
        """
        Get merge events and determine the DC's born through merge events.

        A merge event occurs whenever members of two distinct DC at some time
        point are found together in the same DC one time point later.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_merges`
          - :attr:`~.MajorTrack.community_cby_merges`
          - :attr:`~.MajorTrack.community_dby_merges`

        """
        # no merges possible at the first time step
        self.community_merges = [[]]
        self.community_cby_merges = [[]]
        self.community_dby_merges = [[]]
        for i in range(1, self.length):
            _merges = []
            _newDC_by_merge = []
            _destroyedDC_by_merge = []
            res_pop = self.resident_population(i-1, i)
            for comm, membs in self.comm_members[i].items():
                prev_comms = set(
                    [
                        self.individual_membership[i-1].get(memb, None)
                        for memb in membs
                    ]
                )
                if None in prev_comms:
                    prev_comms.remove(None)
                if len(prev_comms) > 1:
                    _merges.append(
                        (list(prev_comms), [comm])
                    )

                # check if comm is new
                if comm not in prev_comms:
                    # comm is new, check if pure merge or split-merge
                    prev_membs = [
                            pm
                            for prev_comm in prev_comms
                            for pm in self.comm_members[i-1][prev_comm]
                            if pm in res_pop
                            ]
                    if all([pm in membs for pm in prev_membs]):
                        # all members from the merging groups are now here
                        # so this is a pure merge
                        _newDC_by_merge.append(
                                (list(prev_comms), [comm])
                                )
                dest_prev_comms = [
                        pc
                        for pc in prev_comms
                        if pc not in set(self.dcs[i])
                        ]
                for dpc in dest_prev_comms:
                    prev_membs = [
                            pm
                            for pm in self.comm_members[i-1][dpc]
                            if pm in res_pop
                            ]
                    if all([pm in membs for pm in prev_membs]):
                        # prev comm does not exist anymore and remaining membs
                        # are all in comm > merge without split
                        _destroyedDC_by_merge.append(
                                ([dpc], [])
                                )
            self.community_cby_merges.append(_newDC_by_merge)
            self.community_dby_merges.append(_destroyedDC_by_merge)
            self.community_merges.append(_merges)

    def get_community_growths(self,):
        """
        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_growths`
        """
        # birth events are not growth events
        self.community_growths = [[]]
        for i in range(1, self.length):
            _growths = []
            res_pop = self.resident_population(i-1, i)
            for comm, membs in self.comm_members[i].items():
                # make sure it's not another event (birth, split, merge)
                if comm in set(self.dcs[i-1]):
                    new_membs = [memb for memb in membs if memb not in res_pop]
                    if len(new_membs):
                        _growths.append((comm, len(new_membs)))
            self.community_growths.append(_growths)

    def get_community_shrinkages(self,):
        """
        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_shrinkages`
        """
        # birth events are not growth events
        self.community_shrinkages = [[]]
        for i in range(1, self.length):
            _shrinkages = []
            res_pop = self.resident_population(i-1, i)
            for comm, membs in self.comm_members[i-1].items():
                # make sure it's not another event (death, split, merge)
                if comm in set(self.dcs[i]):
                    gone_membs = [mb for mb in membs if mb not in res_pop]
                    if len(gone_membs):
                        _shrinkages.append((comm, len(gone_membs)))
            self.community_shrinkages.append(_shrinkages)

    def get_community_events(self,):
        """
        Compute all dynamic community life-cycle related events.

        """
        self.get_community_births()
        self.get_community_deaths()
        self.get_community_splits()
        self.get_community_merges()
        self.get_community_growths()
        self.get_community_shrinkages()

    def _get_auto_correlation(
        self, community_id, from_idx, to_idx, residents=True
    ):
        r"""
        Get the auto-correlation of members in a dc between two time points.

        The auto-correlation is computed for `community_id` between the time
        points `from_idx` and `to_idx`.
        If the community is not present a both time points None is returned.

        Parameters
        ==========
        community_id: int
          is the identifier of the dynamic cluster.
          .. todo::

            Verify the type

        from_idx: int
          first index of the pairs of slices to compute the auto-correlation
          for.
        to_idx: int
          second index of the pairs of slices to compute the auto-correlation
          for.
        residents: bool (default=True)
          determines if only resident data sources (i.e. that are present in
          both slices) should be considered, or also data sources present in
          only one of the two slices.

        Returns
        =======
        auto_correlation: float
          the autocorrelation given by:

          .. math::

            \frac{|dc_{i} \cap dc_{j}|}{|dc_{i} \cup dc_{j}|_{res}}

          where :math:`i, j` are the indices `from_dix` and `to_idx` and
          :math:`|<selection>|_{res}` is the number of data sources within
          `<selection>` counting all data sources (if `residents=False`) or
          only those present in both slices (`residents=True`).

        """
        res_pop = self.resident_population(from_idx, to_idx)
        if community_id in set(self.dcs[from_idx]) and \
                community_id in set(self.dcs[to_idx]):
            from_comm_membs = self.comm_members[from_idx][community_id]
            to_comm_membs = self.comm_members[to_idx][community_id]
            if residents:
                from_comm_membs = [
                        memb for memb in from_comm_membs
                        if memb in res_pop
                        ]
                to_comm_membs = [
                        memb for memb in to_comm_membs
                        if memb in res_pop
                        ]
            return len(
                    [
                        fm
                        for fm in from_comm_membs
                        if fm in to_comm_membs
                        ]
                    ) / float(len(set(from_comm_membs + to_comm_membs)))
        else:
            return None

    def get_auto_corrs(self, residents=True):
        """
        Get the auto-correlation between any two consecutive slices.

        This method computes for all dynamic clusters the auto-correlation
        between any two consecutive slices, if the dynamic community exists in
        both.
        If residents==True, then only the individuals present in both time
        points are considered.

        Parameters
        ==========
        residents: bool (default=True)
          determines if only resident data sources (i.e. that are present in
          both slices) should be considered, or also data sources present in
          only one of the two slices.

        Returns
        =======
        None: None
          Adds new attributes:

          - :attr:`~.MajorTrack.community_autocorrs`
        """
        self.community_autocorrs = {}
        for idx in range(1, self.length):
            for comm_id in set(self.dcs[idx]):
                auto_corr = self._get_auto_correlation(
                        comm_id, idx - 1, idx, residents
                        )
                if auto_corr is not None:
                    # todo: use defaultdict to avoid handling exception
                    try:
                        self.community_autocorrs[comm_id][idx] = auto_corr
                    except KeyError:
                        self.community_autocorrs[comm_id] = {
                                idx: auto_corr
                                }
        return None

    def get_community_coloring(self, n=None, iterator=None, **kwargs):
        if n is None:
            n = self.comm_nbr
        use_dc = kwargs.get('distinct_colors', None)
        if use_dc is not None:
            # todo: make sure use_dc is an instance of DistinctColors
            assert n <= use_dc.n
            self.color_sequence = use_dc
        elif self.color_sequence is None:
            dc_params = kwargs.get(
                    'dc_params', dict(
                        h_shuffle=False,
                        s_shuffle=False,
                        h_init=0.1,
                        s_init=1.0,
                        v_init=1.0,
                        )
                    )
            self.color_sequence = DistinctColors(
                        n,
                        [0, 1], [0.2, 1.0], [0.9, 1.0],
                        **dc_params
                        )
        else:
            assert self.color_sequence.n >= n
        if iterator is None:
            iterator = range(self.length)
        self.comm_colours = self.color_sequence.get_colors()
        self.sp_next_colour = iter(range(len(self.comm_colours)))
        self.sp_community_colour_idx = {}
        for idx in iterator:
            _grouping = self.clusterings[idx]
            _communities = self.dcs[idx]
            for _i in range(len(_grouping)):
                if _communities[_i] in self.sp_community_colour_idx:
                    g_color_idx = self.sp_community_colour_idx[
                            _communities[_i]]
                else:
                    g_color_idx = next(self.sp_next_colour)
                    self.sp_community_colour_idx[
                            _communities[_i]] = g_color_idx

    def get_alluvialdiagram(
            self, axes, iterator=None, cluster_width=timedelta(days=1),
            *args, **kwargs
            ):
        """
        Takes a matplotlib axes and draws an alluvialdiagram on it.
        `iterator` is the iterator s to draw the clusters for.
        If `iterator is not provided, then the alluvialdiagram will contain all
        the clustrings in the time series.

        Parameters
        ==========
        axes: :obj:`matplotlib.axes.Axes`
          Axes to draw an Alluvial diagram on.
        iterator: iter (default=None)
          An iterator for the indices of the time series to include in the
          alluvial diagram. If not provided then the entire time series is
          used.
        cluster_width: float
          with of the clusters. This should be provided in the same units as
          :attr:`~Majortrack.timepoints`.
        \*args optional parameter:
          Will be forwarded to the :class:`pyalluv.AlluvialPlot` call.
        \**kwargs optional parameter:
          cluster_location: str (default='center')
            either 'center', 'start, 'end' location withing the aggregation
            time window where the cluster should be put.
          cluster_label: str (default=None)
            determine how to label cluster. Possible options are:

            * 'groupsize'
            * 'group_index'

          merged_edgecolor: str (default=None)
            edgecolor of merged clusters.
          merged_facecolor: str (default=None)
            facecolor of merged clusters.
          merged_linewidth: float (default=None)
            linewidth of merged clusters.
          cluster_facecolor: str, dict(dict)
            facecolor of clusters. Either provide a single color or a `dict`
            with indices of the time series as keys, holding a dict with
            cluster_id as key and colours as values.
          cluster_edgecolor: str, dict(dict)
            edgecolor of clusters. Either provide a single color or a `dict`
            with indices of the time series as keys, holding a dict with
            cluster_id as key and colours as values.
          flux_facecolor: str, dict
            either provide a single color, a keyword or a dict.

            Valid keywords are: ``'cluster'``.

            If a dictionary is provided then the `idx` of the time series must
            be the keys with another dict as value holding a dict with a tuple
            as key and a color as value. The tuple's first element must  be a
            group id form time step `idx` and the second a group id k form time
            step `idx`+1
          new_coloring: bool (default=False)
            if a new color sequence should be generated or not.
          distinct_colors: :obj:`colorseq.DistinctColors` (default=None)
            the sequence of distinct colour to use.
          target_clusters: list (default=None)
            list of dynamic cluster id's to display in the alluvial diagram.
            If provided, only the dynamic clusters specified in this list
            will be displayed.
        """
        assert self.clusterings is not None
        # from matplotlib import pyplot as plt
        cluster_label = kwargs.get('cluster_label', None)
        cluster_facecolor = kwargs.get('cluster_facecolor', None)
        cluster_edgecolor = kwargs.get('cluster_edgecolor', None)
        cluster_linewidth = kwargs.get('cluster_linewidth', None)
        cluster_label_margin = kwargs.get('cluster_label_margin', None)
        cluster_location = kwargs.get('cluster_location', 'center')
        flux_facecolor = kwargs.get('flux_facecolor', None)
        flux_edgecolor = kwargs.get('flux_edgecolor', None)
        flux_linewidth = kwargs.get('flux_linewidth', None)
        _dnec = kwargs.get('default_cluster_edgecolor', 'none')
        _dnfc = kwargs.get('default_cluster_facecolor', 'gray')
        # _dnfc = kwargs.get('default_cluster_facecolor', 'gray')
        # _dfec = kwargs.get('default_flux_edgecolor', 'gray')
        # _dffc = kwargs.get('default_flux_facecolor', 'gray')
        _dnlw = kwargs.get('default_cluster_linewidth', 0.0)
        _dflw = kwargs.get('default_flux_linewidth', 0.0)
        if iterator is None:
            _iterator = iter(range(self.length))
        else:
            _iterator = iter(iterator)
        use_community_coloring = True
        # always generate community colors (complete missing dict keys)
        # if isinstance(cluster_facecolor, dict) \
        #         and isinstance(flux_facecolor, dict):
        #     use_community_coloring = False

        # check if coloring is provided
        distinct_colors = kwargs.get('distinct_colors', None)
        if use_community_coloring and kwargs.get('new_coloring', True):
            if self.dcs is None:
                print(
                    "WARNING: The detection of dynamic clusters has not been "
                    "carried out yet. Plotting now will result in default "
                    "colouring of the clusters."
                )
            else:
                if iterator is not None:
                    _all_comms = []
                    for idx in iter(iterator):
                        _all_comms.extend(
                                [com for com in set(self.dcs[idx])]
                                    )
                    n = len(set(_all_comms))
                    self.get_community_coloring(
                            n, iterator=iter(iterator),
                            distinct_colors=distinct_colors
                            )
                else:
                    if self.comm_colours is None:
                        self.get_community_coloring(
                            distinct_colors=distinct_colors
                            )
        _clusters = []
        _fluxes = []
        _clusters_x_pos = []
        # get the colors

        def _c_None(idx, group_id, *args):
            return None

        def _c_comm(idx, group_id):
            if self.dcs is not None and idx < len(self.dcs):
                _c = self.dcs[idx][group_id]
                _c_idx = self.sp_community_colour_idx[_c]
                return self.comm_colours[_c_idx]
            else:
                return _dnec

        def _fc_individ(idx, group_id):
            nfc = cluster_facecolor.get(idx, {}).get(group_id, None)
            return nfc if nfc is not None else _c_comm(idx, group_id)

        def _ec_individ(idx, group_id):
            nec = cluster_edgecolor.get(idx, {}).get(group_id, None)
            return nec if nec is not None else _c_comm(idx, group_id)

        def get_flux_c(flux_fc, flux_ec):
            # for facecolor
            if isinstance(flux_fc, str) and any(
                    cluster in flux_fc for cluster in ['cluster']
                    ):
                def get_flux_fc(idx, g_form, g_to):
                    return flux_fc
            elif isinstance(flux_fc, dict):
                def get_flux_fc(idx, g_from, g_to):
                    ffc = flux_fc.get(idx, {}).get((g_from, g_to), None)
                    if ffc is None:
                        ffc = flux_fc.get(idx, {}).get(g_from, None)
                    return ffc if ffc is not None else _c_comm(idx, g_from)
            else:
                get_flux_fc = _c_None
            # for edgecolor
            if isinstance(flux_ec, str) and any(
                    cluster in flux_ec for cluster in ['cluster']
                    ):
                def get_flux_ec(idx, g_form, g_to):
                    return flux_ec
            elif isinstance(flux_ec, dict):
                def get_flux_ex(idx, g_from, g_to):
                    fec = flux_ec.get(idx, {}).get((g_from, g_to), None)
                    if fec is None:
                        fec = flux_ec.get(idx, {}).get(g_from, None)
                    return fec if fec is not None else _c_comm(idx, g_from)
            else:
                get_flux_ec = _c_None
            return get_flux_fc, get_flux_ec

        def _def_linewidth(*args):
            return _dflw if len(args) == 3 else _dnlw

        def _nlw_individ(idx, group_id):
            return cluster_linewidth.get(idx, {}).get(group_id, _dnlw)

        def _flw_individ(idx, g_from, g_to):
            flw = flux_linewidth.get(idx, {}).get((g_from, g_to), None)
            if flw is None:
                flw = flux_linewidth.get(idx, {}).get(g_from, None)
            return flw if flw is not None else _dflw

        if isinstance(cluster_facecolor, str) or cluster_facecolor is None:
            if cluster_facecolor is None or cluster_facecolor == 'community':
                _get_cluster_facecolor = _c_comm
            else:
                def _get_cluster_facecolor(idx, group_id): return _dnfc
        elif isinstance(cluster_facecolor, dict):
            _get_cluster_facecolor = _fc_individ
        else:
            _get_cluster_facecolor = _c_None

        if isinstance(cluster_edgecolor, str):
            if cluster_edgecolor == 'community':
                _get_cluster_edgecolor = _c_comm
            else:
                def _get_cluster_edgecolor(idx, group_id): return _dnec
        elif isinstance(cluster_edgecolor, dict):
            _get_cluster_edgecolor = _ec_individ
        else:
            _get_cluster_edgecolor = _c_None

        _get_flux_facecolor, _get_flux_edgecolor = get_flux_c(
                flux_facecolor, flux_edgecolor
                )

        if isinstance(cluster_linewidth, dict):
            _get_cluster_linewidth = _nlw_individ
        else:
            _get_cluster_linewidth = _def_linewidth

        if isinstance(flux_linewidth, dict):
            _get_flux_linewidth = _flw_individ
        else:
            _get_flux_linewidth = _def_linewidth

        target_clusters = kwargs.get('target_clusters', False)
        if not target_clusters:
            def include_cluster(idx, group):
                return True
        else:
            def include_cluster(idx, group):
                if any(
                        [
                            indiv in self.comm_members[idx][comm]
                            for indiv in group
                            for comm in target_clusters
                            if comm in self.comm_members[idx]
                            ]
                        ):
                    return True
                elif idx and any(
                        [
                            indiv in self.comm_members[idx-1][comm]
                            for indiv in group
                            for comm in target_clusters
                            if comm in self.comm_members[idx-1]
                            ]
                        ):
                    # check if in the previous snapshot members of this group
                    # were in the target_clusters
                    return True
                else:
                    return False

        merged_edgecolor = kwargs.get('merged_edgecolor', None)
        merged_facecolor = kwargs.get('merged_facecolor', None)
        merged_linewidth = kwargs.get('merged_linewidth', None)
        treat_mergeds = False
        if any([
                    attr is not None
                    for attr in [
                        merged_edgecolor,
                        merged_facecolor,
                        merged_linewidth
                        ]
                    ]) and self.dcs is not None:
            treat_mergeds = True
        for idx in _iterator:
            _grouping = self.clusterings[idx]
            if cluster_location == 'stop':
                # move = self.timepoints[idx][1] - self.timepoints[idx][0]
                move = self.slice_widths[idx]
            elif cluster_location == 'start':
                move = 0
            else:  # this is 'center'
                # move = 0.5 * (
                #     self.timepoints[idx][1] - self.timepoints[idx][0]
                # )
                move = 0.5 * self.slice_widths[idx]
            time_point = self.timepoints[idx] + move
            _clusters_x_pos.append(
                time_point
            )
            tp_clusters = []
            for _i in range(len(_grouping)):
                _group = _grouping[_i]
                if include_cluster(idx, _group):
                    # is it a merged group
                    _facecolor = _get_cluster_facecolor(idx, _i)
                    _edgecolor = _get_cluster_edgecolor(idx, _i)
                    _linewidth = _get_cluster_linewidth(idx, _i)
                    if isinstance(cluster_label_margin, dict):
                        _label_margin = cluster_label_margin[idx][_i]
                    elif isinstance(cluster_label_margin, (tuple, list)):
                        _label_margin = cluster_label_margin
                    else:
                        _label_margin = (0.4 * cluster_width, 0.1)
                    if _edgecolor is None:
                        _edgecolor = _facecolor
                    if treat_mergeds:
                        if idx < len(self.cluster_trace):
                            if self.cluster_trace[idx][_i] < 0:
                                if merged_linewidth is not None:
                                    _linewidth = merged_linewidth
                                if merged_edgecolor is not None:
                                    _edgecolor = merged_edgecolor
                                if merged_facecolor is not None:
                                    _facecolor = merged_facecolor
                    if cluster_label == 'groupsize':
                        c_label = '{0}'.format(len(_group))
                    elif cluster_label == 'group_index':
                        c_label = '{0}'.format(_i)
                    else:
                        c_label = None
                    tp_clusters.append(
                            Cluster(
                                height=len(_group),
                                width=cluster_width,
                                facecolor=_facecolor,
                                edgecolor=_edgecolor,
                                lw=_linewidth,
                                label=c_label,
                                label_margin=_label_margin
                            )
                        )
                else:
                    tp_clusters.append(None)
            _clusters.append(tp_clusters)
        # 2nd round to set the fluxes
        if iterator is None:
            _iterator = iter(range(self.length))
        else:
            _iterator = iter(iterator)
        prev_idx = next(_iterator)
        for idx in _iterator:
            _grouping = self.clusterings[idx]
            _prev_grouping = self.clusterings[prev_idx]
            _x_pos_fluxes = []
            for _i in range(len(_grouping)):
                if include_cluster(idx, _grouping[_i]):
                    for _j in range(len(_prev_grouping)):
                        if include_cluster(idx-1, _prev_grouping[_j]):
                            _intersect = len(
                                    _prev_grouping[_j].intersection(
                                        _grouping[_i]
                                        )
                                    )
                            _linewidth = _get_flux_linewidth(prev_idx, _j, _i)
                            if _intersect:
                                _facecolor = _get_flux_facecolor(
                                        prev_idx, _j, _i)
                                _edgecolor = _get_flux_edgecolor(
                                        prev_idx, _j, _i)
                                if _edgecolor is None:
                                    _edgecolor = _facecolor
                                _x_pos_fluxes.append(
                                        [
                                            Flux(
                                                flux=_intersect,
                                                source_cluster=_clusters[
                                                    prev_idx][_j],
                                                target_cluster=_clusters[
                                                    idx][_i],
                                                facecolor=_facecolor,
                                                edgecolor=_edgecolor,
                                                lw=_linewidth
                                                )
                                            ]
                                        )
            _fluxes.append(
                    _x_pos_fluxes
                    )
            prev_idx = idx
        self.sp_clusters = {
                _clusters_x_pos[idx]: [
                    cluster for cluster in _clusters[idx]
                    if cluster is not None
                    ]
                for idx in range(len(_clusters))
                }
        self.sp_fluxes = {
                _clusters_x_pos[idx]: _fluxes[idx]
                for idx in range(len(_fluxes))
                }

        self.alluvial = AlluvialPlot(
                self.sp_clusters,
                axes,
                *args,
                **kwargs
                )

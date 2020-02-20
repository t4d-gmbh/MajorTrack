#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module documentation goes here
   and here
   and ...
"""


class MajorTrack(object):
    r"""

    Parameters
    ===========

    clusterings: list, dict
      Sequence of clusterings.
      **If provided as a `dict`**:
        keys: float, datetime
          The time points.
        values: list, dict
          The membership list of each clustering indicating to which cluster
          a data source belongs.
          See :obj:`~MajorTrack.memberships` for details.
    \**kwargs optional parameter:
      timepoints: list
        The time points of each clustering.

        Note
        -----
        If `clusterings` if of type `dict` then the keys will be used as time
          points and this optional parameter is ignored, even if provided.

      group_matchup_method: str (default='fraction')
        Set the method to calculate the similarity between two clusters from
        different clusterings. By default the fraction of identical members is
        used as explained in
        `the original article <https://arxiv.org/abs/1912.04261>'_.

    Attributes
    ==========
    group_matchup: list
      Holds for each time point the tracing and mapping sets of all clusters.
      Each element is a `dict` with the keys ``'forward'`` and ``'backward'``.
      Both hold a `dict` indicating for a cluster the best matching cluster
      along with the similarity score of the particular relation in a `tuple`.
      
      Example
      -------
      :: code-block: python

        self.group_matchup[1] = {
            'backward': {0: (0, 1.0), ...},
                         ^   ^  ^
                         |   |  similarity score
                         |   cluster from previous time point
                         cluster from current time point.
            }

    """
    def __init__(self, clusterings, **kwargs):
        assert isinstance(clusterings, (list, dict))
        if isinstance(clusterings, list):
            self.timepoints = kwargs.pop(
                    'timepoints', list(range(len(clusterings)))
                    )
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

        self.group_matchup_method = kwargs.get(
                    'group_matchup_method',
                    'fraction'
                    )

    def get_group_matchup(self, matchup_method=None):
        r"""
        Determine majority relation between neighbouring snapshots.

        Parameters
        ===========
        matchup_method: str (default=None)
          If provided this overwrites `self.group_matchup_method. It determines
          the method to use when calculating similarities between clusters from
          neighbouring snapshots.

        Returns
        =======
        self: :class:`.MajorTrack`
          with new attribute :ref:`group_matchup`.

        ########
        Between each pair of consecutive time points all groups are compared
        and matched (if possible) using `matchup_method`.

        Set:
        ----
        - self.group_matchup: List holding for each time point a dictionary
            with 'backward'/'forward' matchups. A matchup is a dict indicating
            for each group (id) the best match. The best match is given by a
            tuple with group id and similarity score.
            E.g.: self.group_matchup[1] = {
                'backward': {0: (0, 1.0), ...},
                'forward': {0: (1, 0.7), ...}
                }
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
                        for _group_id in range(len(self.groupings[0]))
                        }
                    }
                )
        self.group_similarities.append(
                {
                    'backward': {
                        _group_id: None
                        for _group_id in range(len(self.groupings[0]))
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
                for _group_id in range(len(self.groupings[-1]))
                }
        self.group_similarities[-1]['forward'] = {
                _group_id: None
                for _group_id in range(len(self.groupings[-1]))
                }

    def get_span(self, idx, span_set, get_indivs=True):
        r"""
        Create the tracer tree.

        Parameters
        ===========

        ####
        Get the span (time forward)
        """
        span_tree = {}
        if isinstance(span_set, int):
            span_tree[idx] = [self.groupings[idx][span_set]]
        elif isinstance(span_set, str):
            span_tree[idx] = filter(
                    lambda g: span_set in g,
                    self.groupings[idx]
                    )
        else:
            span_tree[idx] = [span_set]
        current_set = set.union(*span_tree[idx])
        for _idx in range(idx + 1, self.length):
            next_groupings = self.groupings[_idx]
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
                        i for i in range(len(self.groupings[_idx]))
                        if any([
                            el in _span_set
                            for el in self.groupings[_idx][i]
                            ])
                        ]
            return span_tree_idxs

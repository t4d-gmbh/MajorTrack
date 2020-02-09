from __future__ import absolute_import, \
        print_function, unicode_literals, division


def node_l_index(node_id, *groupings, **kwargs):
    """
    Return the node L-index for a node and a set of groupings.

    .. todo::

        Add option to only consider nodes that ever were in the same group
        as the focal node.



    Parameters
    ----------
    node_id: int
      Identifier of the focal node.
    groupings: list
      Each grouping must be a list of lists where each element is a node.
      So a grouping must be of the form:

      .. code-block:: python

            [[n0, node_id, n2], [n3, n4, ...], ...]
            # <    group 1  >    < group 2 > ...

      .. note::

          Each grouping must contain the same set of nodes.
    """
    nbr_groupings = kwargs.get('nbr_groupings', len(groupings))
    nbr_nodes = kwargs.get(
            'nbr_nodes',
            sum([len(_g) for _g in groupings[0]])
            )
    # TODO: implement the get all nodes in a single statement
    _others = kwargs.get('all_nodes', [])
    # get set of other nodes
    _others = []
    _to_del = [_others.extend(_group) for _group in groupings[0]]
    del _to_del
    _others.remove(node_id)
    # get indexes of group containing node_id in each grouping
    _foc_group_ids = []
    for grouping in groupings:
        _foc_group_ids.append(
                next(
                    idx
                    for idx in range(len(grouping))
                    if node_id in grouping[idx]
                    )
                )
    # get the per node 'energy'
    per_node_count = [
            # for each _grouping:
            [
                1 if _o_node in groupings[i][_foc_group_ids[i]] else - 1
                for _o_node in _others
                ]
            for i in range(len(groupings))
    ]
    # sum over the groupings for each node
    # and normalize by the number of groupings
    # absolute value
    per_node_e = [
            abs(_per_node_sum / nbr_groupings)
            for _per_node_sum in map(sum, zip(*per_node_count))
            ]
    # sum and normalize by the number of other nodes
    return 1 / (nbr_nodes - 1) * sum(
            per_node_e
            )


def l_index(*groupings):
    """
    .. todo::

        This is not implemented.

    """
    # nbr_groupings = len(groupings)
    # nbr_nodes = sum([len(_g) for _g in groupings[0]])
    # get set of other nodes
    _others = []
    _to_del = [_others.extend(_group) for _group in groupings[0]]
    del _to_del


def group_similarity_fraction(from_group, to_group):
    """
    Fraction of nodes present in one group that come from the another group.

    Parameters
    ==========

    from_group: set
      group from which nodes might come.
    to_group: set
      group in which the nodes are.
    """
    if len(to_group):
        return len(to_group.intersection(from_group))/len(to_group)
    else:
        return 0.0


def group_similarity_jaccard(group_1, group_2, N=None):
    """
    Compute the similarity, i.e. the fraction of agreeing nodes, between two
    groups of different groupings.

    The similarity is given by 1 - the fraction of nodes belonging to the
    symmetric difference (or disjunctive union) between the two groups.
    Or, said differently, the fraction of nodes that agree on their respective
    status towards both provided groups. The status of a node towards a group
    can either be member or non-member of that group. A node with identical
    status towards both groups contributes with a summand of 1/N to the
    similarity of those two groups.

    If the total size of the networks, N, is provided, then also nodes
    contribute to the similarity of the groups which are non-members for both
    of the groups. Without the total size of the network, `N=None`, the
    groups similarity only considers nodes present in both groups as
    contributing to the similarity and thus is equivalent to the Jaccard index.

    Parameter
    ---------
    group_1: set
      members of a group
    group_2: set
      members of another group
    N: int (default=None)
      Number of nodes in the network. If not provided, then the Jaccard index
      is computed.
    """
    if N is None:
        N = len(group_1.union(group_2))
    if N:
        return (N - len(group_1.symmetric_difference(group_2))) / N
    else:
        return 0.0

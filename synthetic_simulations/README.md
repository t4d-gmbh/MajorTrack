# Generation of community-structured time-series synthetic data

The code provided generates a set of `R` matrices describing the interaction of `N` agents throughout `T` time steps. Agents randomly interact in evenly distributed time intervals `\DeltaT_k` according to specified probability matrices `P_k`. These probability matrices have a community structure, meaning that agents are divided in groups, and agents within each group have a higher probability of interaction with members of its own group than with members of other groups.


### Scripts

The code has the following scripts and functions:

* `generate_partitions`: Given a set of `N` agents and `K` groups of sizes specified in a vector of integers such that `sum(subset_sizes) = N`, it generates a vector randomly associating to each agent a group membership id. It returns the vector `assignments`.

* `modify_partitions`: Modifies a vector of assignment according to rules specified in a list termed `operation`. The operation list should be specified as follows depending on the intended operation:
   * Merging groups: It merges the groups specified into a single group. The list
     should contain:
       * `operation$type = "merge"`. A string specifying the type of operation
       * `operation$groups`: A vector of integers specifing the groups to be merged into the
          one specified as first element of the vector.
   * Splitting a group: It split a groups into two or more subgroups of specified size. The list elements should be
       * `operation$type = "split"`. A string specifying the type of operation
       * `operation$group`: An integer specifying the group to split.
       * `òperation$sizes`: A vector of integers specifying the sizes of the new groups, being the first element the remaining size of the original group.
   * Drift between groups: It moves a fixed number of agents between the group specified and in the order in which they appear. For instance it could move to elements from group 1 to 2, another 2 (different ones)  from 2 to 3, and so on.
      * `operation$type = "drift"`. 
      * `operation$groups`: An ordered vector of integers specifying the groups transfering their elements.
      * `operation$size`: An integer indicating the number of elements transferred between groups.

* `simulate_interactions`:This function randomly draws time_steps matrices of interaction between `N` elements according to the the probabilities specified by probability matrices contained in the list `probability\_matrices`. Therefore, the total number of matrices generated will be `time\_steps * length(probabiliy\_matrices)`.

* `example_use`: A script with a minimal working example describing the use of the above scripts.

### Simulations and output data

The simulations presented in the paper were coded in the script `main`. It requires three bits of code in which the operations made on the groups were specified in three increasingly complex problems, in the files `problem_merge.R`, `problem_merge-split.R` and `problem_merge-split-drift.R`. Each of them contains the operations indicated, and they share some of the starting operations. Since the seed was set, it is possible to observe the effects of the different operations.

We considered three increasingly challenging operations, in which the probability of observing one interaction at a given time step was fixed to `P_between = 0.05`, while the probability of observing one interaction within groups was decreasing: `P_within = 0.8, 0.5, 0.2`. Therefore, a total of 6 scenarios (3 problems and 3 combinations of probabilities) were analysed.

Each run of code generate a folder named with a label `$label = $problem\_$parameters` describing the problem and parameters considered, containing:

* Several heatmaps plotting the probability matrices.
* One randomly-selected interaction matrix for each probability matrix. 
* `matrices_interaction_$label.RDS`: This is a list with the set of interaction matrices generated. This should be the input for the community detection algorithm.
* `assignments_$label.RDS`: This is a list of vectors with the assignments from which the probability matrices (and, in turn, the interaction matrices) were generated. It should be considered the ground truth against which testing the community detection output.



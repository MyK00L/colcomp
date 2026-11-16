# Some other things
## Preprocessing
### Extended rule 4:
When considering an edge for deletion, instead of looking at a single component, you can look at both ends.
If the minimum cut that isolates those components is $\le$ than the sum of the cuts between nodes of the same
colors on both of the components, you can remove the edge.

Furthermore, you can consider any amount of components:
let's say the 2 components on the sides of the edge are both only connected to another component,
with each of those edges having both nodes share a color on the other component,
and that component is not connected to anything, you could add its edge-connectivity to the two cuts.
You could also take into consideration many components, and this would become a small optimization problem
that includes the Rule 4 as explained in the report
(represent each component as a node with positive weight equal to the connectivity, or mincut for the two components with the considered edge, and make edges between them with negative weight equal to the number of edges between them in the original graph,
find a subset with positive weight, considering all edges even to the rest of the graph).
This would probably not be any better on random graphs though.

Also, considering multiple edges between two components, you could generate some cuts that might be useful in the pricing subproblem.

### Preprocessing during the BnP
Right now the preprocessing (like the heuristics) works on a `problem`, if it is made to work on an `mgraph`,
it could be used inside the _BnP_ too.
Also `local_search` could be used to improve primal bounds.

## Cuts
As mentioned in the report, cuts (eg clique and odd-cycle inequalities) could be added to the master problem.

## Errata
In the report it is mentioned that the large scale neighbourhood search did not improve a lot in smaller graphs.
This is because it was not implemented properly at the time and it would not consider paths ending by adding a node
to a component that did not contain its color.


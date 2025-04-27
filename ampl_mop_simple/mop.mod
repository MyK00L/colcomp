param n;
param m;

set V := 1..n;
set E within {u in V, v in V: u<v};
set P := 1..n; # partitions

param col{V}; # color of each node

var x{P cross V} binary; # partition of nodes
var c{P cross E}; # auxilliary variable c[p][i][j] = x[p][i]&x[p][j]

subject to c_and_u{p in P, (u,v) in E}:
	c[p,u,v] <= x[p,u];
subject to c_and_v{p in P, (u,v) in E}:
	c[p,u,v] <= x[p,v];
subject to col_diff{p in P, (u,v) in {u in V, v in V: u<v}: col[u]=col[v]}:
	x[p,u]+x[p,v]<=1;
subject to partition{u in V}:
	sum {p in P} x[p,u]=1;

maximize z: sum{p in P, (u,v) in E}c[p,u,v];


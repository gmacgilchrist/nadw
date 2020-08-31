function D=ariane_D(distchar,edges,S)

[D.n,D.edges,D.bins]=histcounts(distchar(S),edges);
D.S=S;
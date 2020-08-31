function D=ariane_D2(distchar1,distchar2,edges1,edges2,S)

[D.n,D.edges1,D.edges2,D.bins1,D.bins2]=histcounts2(distchar1(S),distchar2(S),edges1,edges2);
D.S=S;
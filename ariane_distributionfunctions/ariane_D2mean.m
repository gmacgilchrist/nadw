function Dmean = ariane_D2mean(meanchar, D)

C = meanchar(D.S);

Dmean = nan(length(D.edges1)-1,length(D.edges2)-1);

BINS = unique([D.bins1 D.bins2],'rows');
if issame(BINS(1,:),[0 0]); BINS(1,:)=[]; end
for b = 1:size(BINS,1);
    i=BINS(b,1); j=BINS(b,2);
    ind = D.bins1==i & D.bins2==j;
    if sum(ind)>0;
        Dmean(i,j)=mean(C(ind));
    end
end
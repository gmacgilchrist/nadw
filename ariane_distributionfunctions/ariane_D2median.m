function Dmedian = ariane_D2median(medianchar, D)

C = medianchar(D.S);

Dmedian = nan(length(D.edges1)-1,length(D.edges2)-1);

BINS = unique([D.bins1 D.bins2],'rows');
if issame(BINS(1,:),[0 0]); BINS(1,:)=[]; end
for b = 1:size(BINS,1);
    i=BINS(b,1); j=BINS(b,2);
    ind = D.bins1==i & D.bins2==j;
    if sum(ind)>0;
        Dmedian(i,j)=median(C(ind));
    end
end
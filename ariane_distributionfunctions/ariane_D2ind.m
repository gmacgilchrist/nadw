function Dind = ariane_D2ind(D)

Dind = cell(length(D.edges1)-1,length(D.edges2)-1);

BINS = unique([D.bins1 D.bins2],'rows');
index = 1:length(D.bins1);
if issame(BINS(1,:),[0 0]); BINS(1,:)=[]; end
for b = 1:size(BINS,1);
    i=BINS(b,1); j=BINS(b,2);
    ind = D.bins1==i & D.bins2==j;
    if sum(ind)>0;
        Dind{i,j}=index(ind);
    end
end
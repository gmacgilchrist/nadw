function [Dsum,Dind] = ariane_Dsum(sumchar, D)

C = sumchar(D.S);

Dsum = nan(length(D.edges)-1,1);
Dind = cell(length(D.edges)-1,1);
index = 1:length(D.bins);

BINS = unique(D.bins);
if BINS(1)==0; BINS(1)=[]; end
for b = 1:size(BINS,1);
    i=BINS(b,1);
    ind = D.bins==i;
    Dind{b}=index(ind);
    Dsum(b)=sum(C(ind));
end
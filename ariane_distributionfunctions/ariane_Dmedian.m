function Dmedian = ariane_Dmedian(medianchar, D)

C = medianchar(D.S);

Dmedian = nan(length(D.edges)-1,1);
BINS = unique(D.bins);
if BINS(1)==0; BINS(1)=[]; end
for b = 1:size(BINS,1);
    i=BINS(b,1);
    ind = D.bins==i;
    Dmedian(b)=median(C(ind));
end
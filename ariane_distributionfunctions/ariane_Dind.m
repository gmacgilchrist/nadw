function Dind=ariane_Dind(D)
% function Dind = ariane_Dindex(opchars,D,operation)
%
% Return the indices of the particles in each histogram bin of D
%
% G.A. MacGilchrist (27/03/18) gmacgilchrist@gmail.com

Dind = cell(length(D.edges)-1,1);
BINS = unique(D.bins);
index = 1:length(D.bins);
if BINS(1)==0; BINS(1)=[]; end
for b = 1:size(BINS,1);
    i=BINS(b,1);
    ind = D.bins==i;
    Dind{b}=index(ind);
end
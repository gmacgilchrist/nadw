function [Dind,Dsum,Dmean,Dmedian,Dstd]=ariane_Dop(opchar,D)
% function Dop = ariane_Dop(opchars,D,operation)
%
% Perform operations (sum, average, indexing) on the characteristics in
% each bin of a distribution, defined by the distribution structure array D.
%
% INPUT
%   opchar      ariane output of characteristic on which to operate
%   D           structure array of distribution
%
% OUTPUT
%   Dind        A n x m cell (where n and m are the dimensions of the
%               distribution) with the indices of the particles that fall
%               in each bin.
%   Dsum        A n x m matrix with the sum of opchar for particles in each
%               bin
%   Dmean, Dmedian  As for Dsum but taking the mean and median
%

% G.A. MacGilchrist (27/03/18) gmacgilchrist@gmail.com
disp('ariane_Dop: Indexing, summing and averaging characteristics')
disp('            within each bin of the distribution D');
if isvector(D.n);
    % A 1-dimensional distribution
    C = opchar(D.S);
    
    index = 1:length(D.bins);
    Dind = cell(length(D.edges)-1,1);
    Dsum = nan(length(D.edges)-1,1);
    Dmean = nan(length(D.edges)-1,1);
    Dmedian = nan(length(D.edges)-1,1);
    
    BINS = unique(D.bins);
    if BINS(1)==0; BINS(1)=[]; end
    n=floor(log(abs(size(BINS,1)))./log(10)); % Order of magnitude
    for b = 1:size(BINS,1);
        i=BINS(b,1);
        ind = D.bins==i;
        Dind{i}=index(ind);
        Dsum(i)=sum(C(ind));
        Dmean(i)=mean(C(ind));
	Dstd(i)=std(C(ind));
        Dmedian(i)=median(C(ind));
        
	if mod(b,10^(n-1))==0;
                disp([num2str((b/size(BINS,1))*100) '% complete']);
        end
    end
    
else
    % A 2-dimensional distribution
    C = opchar(D.S);
    
    index = 1:length(D.bins1);
    Dind = cell(length(D.edges1)-1,length(D.edges2)-1);
    Dsum = nan(length(D.edges1)-1,length(D.edges2)-1);
    Dmean = nan(length(D.edges1)-1,length(D.edges2)-1);
    Dmedian = nan(length(D.edges1)-1,length(D.edges2)-1);
    
    BINS = unique([D.bins1 D.bins2],'rows');
    if issame(BINS(1,:),[0 0]); BINS(1,:)=[]; end
    n=floor(log(abs(size(BINS,1)))./log(10));
    for b = 1:size(BINS,1);
        i=BINS(b,1); j=BINS(b,2);
        ind = D.bins1==i & D.bins2==j;
        if sum(ind)>0;
            Dind{i,j}=index(ind);
            Dsum(i,j)=sum(C(ind));
            Dmean(i,j)=mean(C(ind));
	    Dstd(i,j)=std(C(ind));
            Dmedian(i,j)=median(C(ind));
        end

	if mod(b,10^(n-1))==0;
		disp([num2str((b/size(BINS,1))*100) '% complete']);
	end
    end    
end

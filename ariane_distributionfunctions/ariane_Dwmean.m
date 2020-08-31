function [Dwmean,Dwstd]=ariane_Dwmean(opchar,volume,D)
% function Dop = ariane_Dwmean(opchars,volume,D)
%
% Calculate weighted means and standard deviations on the characteristics in
% each bin of a distribution, defined by the distribution structure array D.
%
% INPUT
%   opchar      ariane output of characteristic on which to operate
%   volume      ariane input with initial particle volume
%   D           structure array of distribution
%
% OUTPUT
%   Dwmean      A n x m matrix with the volume weighted mean of opchar for
%               particles in each bin of F
%   Dwstd       As for Dwmean but weighted standard deviation.
%

% G.A. MacGilchrist (25/8/18) gmacgilchrist@gmail.com
disp('ariane_Dop: Volume-weighted average of characteristics')
disp('            within each bin of the distribution D');
if isvector(D.n);
    % A 1-dimensional distribution
    C = opchar(D.S);
    
    index = 1:length(D.bins);
    Dwmean = nan(length(D.edges)-1,1);
    Dwstd = nan(length(D.edges)-1,1);
    
    BINS = unique(D.bins);
    if BINS(1)==0; BINS(1)=[]; end
    n=floor(log(abs(size(BINS,1)))./log(10)); % Order of magnitude
    for b = 1:size(BINS,1);
        i=BINS(b,1);
        ind = D.bins==i;
        Dwmean(i)=sum(volume(ind).*C(ind))./sum(volume(ind));
        Dwstd(i)=sqrt(sum(volume(ind).*(C(ind)-Dwmean(i)).^2)./sum(volume(ind)));
        
	if mod(b,10^(n-1))==0;
                disp([num2str((b/size(BINS,1))*100) '% complete']);
        end
    end
    
else
    % A 2-dimensional distribution
    C = opchar(D.S);
    
    index = 1:length(D.bins1);
    Dwmean = nan(length(D.edges1)-1,length(D.edges2)-1);
    Dwstd = nan(length(D.edges1)-1,length(D.edges2)-1);
    
    BINS = unique([D.bins1 D.bins2],'rows');
    if issame(BINS(1,:),[0 0]); BINS(1,:)=[]; end
    n=floor(log(abs(size(BINS,1)))./log(10));
    for b = 1:size(BINS,1);
        i=BINS(b,1); j=BINS(b,2);
        ind = D.bins1==i & D.bins2==j;
        if sum(ind)>0;
	    Dwmean(i,j)=sum(volume(ind).*C(ind))./sum(volume(ind));
            Dwstd(i,j)=sqrt(sum(volume(ind).*(C(ind)-Dwmean(i)).^2)./sum(volume(ind)));
        end

	if mod(b,10^(n-1))==0;
		disp([num2str((b/size(BINS,1))*100) '% complete']);
	end
    end    
end

function [m,s] = ariane_meanbin(Var,Xb,Xe,Yb,Ye)

switch nargin
    case{3}
        m = nan(length(Xe)-1,1);
        s = nan(length(Xe)-1,1);
        BINS = unique(Xb);
        if BINS(1)==0; BINS(1)=[]; end
        for b = 1:size(BINS,1);
            if mod(b,1000)==0; disp([num2str(b) ' of ' num2str(length(BINS))]); end
            i = BINS(b);
            ind = Xb==i;
            m(i)=mean(Var(ind));
            s(i)=std(Var(ind));
        end
        
    case{5}
        BINS = unique([Xb(:) Yb(:)],'rows');
        if issame(BINS(1,:),[0 0]); BINS(1,:)=[]; end
        m = nan(length(Xe)-1,length(Ye)-1);
        s = nan(length(Xe)-1,length(Ye)-1);
        for b = 1:size(BINS,1);
            if b==1 || mod(b,1000)==0; disp([num2str(b) ' of ' num2str(length(BINS))]); end
            i=BINS(b,1); j=BINS(b,2);
            ind = Xb==i & Yb==j;
            if sum(ind(:))>0;
            m(i,j)=mean(Var(ind));
            s(i,j)=std(Var(ind));
            end
        end
        
end
        
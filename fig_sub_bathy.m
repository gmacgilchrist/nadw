close all
universal
savefigs='y';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VENTILATION AS FUNCTION OF BATHYMETRY AND DENSITY
% CALCULATIONS
if donebathycalcs=='n';
	% BATHYMETRIC CONTOUR OF SUBDUCTION
	
	if exist('matfiles/D2_fXfY_sum_volume.mat')~=2;
        	% Calculate subduction volume for each grid square
       		S = ariane_S({init_t,final_section,final_age},{4217,7,[31536000 Inf]});

        	% Distribution at subduction location
        	minX = floor(min(final_x(S))); maxX = ceil(max(final_x(S)));
        	minY = floor(min(final_y(S))); maxY = ceil(max(final_y(S)));
        	edgesX = minX:1:maxX;
        	edgesY = minY:1:maxY;
        	D = ariane_D2(final_x,final_y,edgesX,edgesY,S);
        	[Dind,Dsum,~,~] = ariane_Dop(init_volume,D);
        	Dsum(Dsum==0)=NaN;
        	save('matfiles/D2_fXfY_sum_volume.mat','D','Dsum','Dind');
        else
        	load('matfiles/D2_fXfY_sum_volume.mat');
        end

	% Assign model bathymetry level to each particle (in subset S)
	final_bathy = nan(sum(S),1);
	[nx,ny]=size(Dind);
	for i=1:nx;
	    for j=1:ny;
	        if ~isempty(Dind{i,j})
	            final_bathy(Dind{i,j}(:))=mbathy(TrangeX(i),TrangeY(j));
	        end
	    end
	end
	% Calculate distribution with respect to model bathymetry level
	Sb = S(S);
	Db = ariane_D(final_bathy,0.5:1:75.5,Sb);
	
	[Dbind,Dbsum,~,~] = ariane_Dop(init_volume(S),Db);
	
	Dbsum(isnan(Dbsum))=0;
	Dbsumcum = cumsum(Dbsum)./sum(Dbsum);
	
	resetcalcs
	donebathycalcs='y';
end

% FIGURE
f1 = figure
subplot(211)
bar(Dbsumcum,1,'facecolor',0.6*[1 1 1],'edgecolor','none');
set(gca,'xlim',[0 75],'ytick',0:0.25:1);
hold on;
% Plot vertical lines for contours in map
ylims = get(gca,'ylim');
plot([bathys; bathys]+0.5,repmat(ylim',[1 length(bathys)]),'-','color',0.8*[1 1 1])
% Sum percentages between contours
fs = 8;
for c = 1:length(bathys);
    if c==1;
        text(bathys(c)-5+0.5,0.10,...
            [num2str(round(abs(diff([Dbsumcum(bathys(c)) 0]))*100)) '\%'],...
            'fontsize',fs)
        text(bathys(c)+0.5,0.90,num2str(round(gdept(bathys(c)),-1)),'fontsize',fs,'horizontalalignment','c')
    else
        text(mean([bathys(c) bathys(c-1)])+0.5,0.10,...
            [num2str(round(abs(diff([Dbsumcum(bathys(c)) Dbsumcum(bathys(c-1))]))*100)) '\%'],...
            'fontsize',fs, 'color','w','horizontalalignment','c')
        text(bathys(c)+0.5,0.90,num2str(round(gdept(bathys(c)),-1)),'fontsize',fs,'horizontalalignment','c')
    end
end
text(25,0.90,'Depth (m) :')
xlabel('Model z-level'); ylabel('Fraction');
title('Bathymetric model level over which subduction occurs','fontsize',titlefs)
if strcmp(savefigs,'y')
    %export_fig([savedir 'fig_subduction_bathy-dens'],'-png',saveres);
%	saveas(f1,[savedir 'fig_subduction_bathy'],'png');
	print(f1, [savedir 'fig_subduction_bathy'], '-dpng',saveres);
end

% DENSITY DISTRIBUTION IN EACH BATHYMETRY BAND
cm = cbrewer('qual','Set1',9);

f2 = figure
rr = 27.7:0.01:28;
% Initial density
Dbid = ariane_D2(final_bathy,init_dens(S),[-Inf bathys Inf],rr,Sb);
[~,Dbidsum,~,~] = ariane_Dop(init_volume(S),Dbid);
Dbfd = ariane_D2(final_bathy,final_dens(S),[-Inf bathys Inf],rr,Sb);
[~,Dbfdsum,~,~] = ariane_Dop(init_volume(S),Dbfd);
subplot(122); hold on
h=cell(3,1);
s=cell(3,1);
for b=2:4;
    pid = Dbidsum(b,:);
    h{b-1}=plot(pid*1e-14,centre(rr),'-','color',cm(b,:),'linewidth',2);
    pfd = Dbfdsum(b,:);
    plot(pfd*1e-14,centre(rr),'--','color',cm(b,:));
    s{b-1}=[num2str(round(gdept(bathys(b-1)),-1)) ' to ' num2str(round(gdept(bathys(b)),-1)) ' m'];
end
xlims = [0 9];
set(gca,'box','on','xlim',xlims,'ydir','r')
ylabel('Neutral density, $\gamma^n$'); xlabel('Volume ($10^{14}m^3$)');
title('Density distribution','fontsize',titlefs)

plot([xlims(1)+2 xlims(1)+3],[27.71 27.71],'k--'); text(xlims(1)+3.5,27.71,'Density at subduction');
plot([xlims(1)+2 xlims(1)+3],[27.73 27.73],'k-','linewidth',2); text(xlims(1)+3.5,27.73,'Density in Sep. 2015');

l = legend([h{1} h{2} h{3}],{s{1},s{2},s{3}});
set(l,'box','off','location','southeast','interpreter','latex')

if strcmp(savefigs,'y')
    %export_fig([savedir 'fig_subduction_bathy-dens'],'-png',saveres);
%	saveas(f2,[savedir 'fig_subduction_dens'],'png');
	print(f2, [savedir 'fig_subduction_dens'], '-dpng',saveres);
end


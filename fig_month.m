close all
universal
savefigs='y';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MONTH OF VENTILATION

%% CALCULATIONS
if donemonthlycalcs=='n'
	disp('PERFORMING MONTHLY CALCULATIONS')
	% Lagrangian
	S = ariane_S({init_t,final_section,final_age},{4217,7,[31536000 Inf]});
	edges = 0+1/(fm*2):1/fm:12+1/(fm*2);
	D = ariane_D(final_month,edges,S);
	[~,Dsum,~,~] = ariane_Dop(init_volume,D);
	% Eulerian
	lookfile = 'matfiles/month_eulerian.mat';
	if exist(lookfile)~=2;
	years = 1975:2015;
	fract = nan(1442,1021,length(years));
	for y = 1:length(years);
	    disp(y)
	    files = dir([modeldir num2str(years(y)) '/ORCA025.L75-GJM189_*_gridT.nc']);
	    somxl010 = nan(1442,1021,length(files));
	    for f = 1:length(files)
	        somxl010(:,:,f) = ncread([modeldir num2str(years(y)) '/' files(f).name ],'somxl010');
	    end
	    [~,time] = max(somxl010,[],3);
	    time(time==1)=NaN;
	    fract(:,:,y)=time/length(files);
	end
		fractNA = fract(region_limits(1,1):region_limits(1,2),region_limits(2,1)+300:region_limits(2,2),:);
		%  Load area file
		area = ncread([rootdir 'grids_ncmod/' gridf ],'e1t').*ncread([rootdir 'grids_ncmod/' gridf ],'e2t');
		areaNA = area(region_limits(1,1):region_limits(1,2),region_limits(2,1)+300:region_limits(2,2));
		areaNA = repmat(areaNA,[1,1,size(fractNA,3)]);
		edges_eul = edges/12;
		S_eul = true(length(fractNA(:)),1);
		D_eul = ariane_D(fractNA(:),edges_eul,S_eul);
		[~,Dsum_eul,~,~] = ariane_Dop(areaNA(:),D_eul);
		save(lookfile,'D_eul','Dsum_eul');
	else
		load(lookfile);
	end
	resetcalcs
	donemonthlycalcs='y'

end

%% FIGURE
close all
% Figure properties
lmap = 0.1; bmap = 0.5; wmap = 0.5; hmap = 0.4;
gapv = 0.06;
hhist = 0.1; lhist = lmap;  bhist = bmap-gapv-hhist; whist = wmap;
hcbar = 0.01; lcbar = lmap; bcbar = bhist-gapv/4-hcbar; wcbar = wmap;
posmap = [lmap bmap wmap hmap];
poshist = [lhist bhist whist hhist];
poscbar = [lcbar bcbar wcbar hcbar];
clims = [0 12]; nc = 12; cm = cbrewer('qual','Paired',nc);
cm = [unventcolor; cm];

f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'pos',[pos(1) pos(2) 8.3 11.7]);

% Map
% Grid initial locations
ind = sub2ind(size(tmask),ceil(init_x),ceil(init_y),floor(init_z));
month = nan(size(tmask));
month(ind(S)) = final_month(S);
month_med = nanmedian(month,3);
month_med(isnan(month_med) & seeded~=0) = 0;

axes('position',posmap);
worldmap(latlimB,lonlimB);
pcolorm(double(gphit),double(glamt),ceil(month_med)); shading flat
framem('FlineWidth',2,'FEdgeColor','black');
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none');
colormap(cm); caxis(clims);

title('(a) Month of subduction in NADW (depth-median)','fontweight','normal','interpreter','none','fontsize',titlefs)

% Histogram
% Lagrangian
cmc = cbrewer('qual','Paired',nc);
axes('position',poshist);
for mt = 1:length(Dsum);
    fc = cmc(ceil(mt/fm),:);
    bar(edges(mt),Dsum(mt)*1e-15,1/fm,'facecolor',fc,'edgecolor','none');
    hold on
end
set(gca,'xticklabels',[]);
box(gca,'off')
xlim([0 12])
ylabel('Volume ($10^{15}\, m^3$)');
% Eulerian
ax_eul = axes('position',poshist);
%for mt = 1:12;
%    fc = cmc(mt,:);
%    if mt==12
%        ind = (mt-1)*fm+1:(mt-1)*fm+fm;
%    else
%        ind = (mt-1)*fm+1:(mt-1)*fm+fm+1;
%    end
%    plot(edges_eul(ind),Dsum(ind)/sum(Dsum),'color','k','linestyle','--')
%    hold on
%end
edges_eul = D_eul.edges;
plot(centre(edges_eul),smooth(Dsum_eul/sum(Dsum_eul)),'linestyle','--','color','r');
set(ax_eul,'color','none','ycolor','r','yaxislocation','right','xticklabels',[])
box(ax_eul,'off')
xlim([0 1])
ylabel('Fraction','interpreter','none')
title({'(b) Volume of NADW subducted in each month','and \color{red}month of deepest mixed layers'},'fontsize',titlefs,'fontweight','normal','interpreter','tex')

% Colour bar
yc = 0.5:12.5;
xc = 0:2;
[~,zc] = meshgrid(xc,yc);
ax = axes('position',poscbar);
contourf(yc,xc,zc',yc,'linestyle','none'); caxis(clims); colormap(ax,cbrewer('qual','Paired',12));
set(gca,'ytick',[],'xtick',1:12,'ticklength',[0 0],...
    'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
xlabel('Month of subduction','interpreter','none')

if strcmp(savefigs,'y')
   % export_fig([savedir 'fig_month'],'-png',saveres);
%	saveas(f,[savedir 'fig_month'],'png');
	print(f, [savedir 'fig_month'], '-dpng',saveres);
end


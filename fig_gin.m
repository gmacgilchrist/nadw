close all
universal
savefigs='y';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GIN SEAS RE-ENTRAINMENT
% WILL RE-RUN EXPERIMENT IN GLOBAL ORCA 5D, SEEDED IN SEPTEMBER

if donegincalcs=='n';

	% Load data from GIN seas forwards-in-time experiment
	
	final_x_g = ncread(fileout_g,'final_x');
	final_y_g = ncread(fileout_g,'final_y');
	init_x_g = ncread(fileout_g,'init_x');
	init_y_g = ncread(fileout_g,'init_y');
	init_volume_g = ncread(filein_g,'init_volume');
	init_dens_g = ncread(filein_g,'init_dens');
	final_dens_g = ncread(fileout_g,'final_dens');
	final_section_g = ncread(fileout_g,'final_section');
	final_lat_g = ncread(fileout_g,'final_lat');
	final_age_g = ncread(fileout_g,'final_age');
	
	glamt_g = ncread([rootdir 'grids_ncmod/' gridf_g ],'glamt');
	gphit_g = ncread([rootdir 'grids_ncmod/' gridf_g ],'gphit');
	[nx_g,ny_g]=size(glamt_g);

	% Calculations
	% Re-entrainment
	Smap = ariane_S({final_section_g,final_age_g,init_x_g,init_y_g},{7,[31536000/4 Inf],[300 Inf],[440 Inf]});
	edgesX = 0:nx_g;
	edgesY = 0:ny_g;
	D = ariane_D2(final_x_g,final_y_g,edgesX,edgesY,Smap);
	[~,Dsum_map,~,~] = ariane_Dop(init_volume_g,D);
	% Fraction of volume at each latitude
	S = ariane_S({final_age_g},{[31536000/4 Inf]});
	edges_lat = latlimB(1):2:latlimB(2);
	D = ariane_D(final_lat_g,edges_lat,S);
	[~,Dsum_lat,~,~] = ariane_Dop(init_volume_g,D);
	Dsum_frac = Dsum_lat/sum(init_volume_g(S)); % Build fraction functionality into ariane_Dop
	Dsum_frac(isnan(Dsum_frac))=0;
	% Fraction of re-entrainment at each latitude
        S_reent = ariane_S({final_age_g,final_section_g},{[31536000/4 Inf],7});
        edges_lat = latlimB(1):2:latlimB(2);
        D = ariane_D(final_lat_g,edges_lat,S_reent);
        [~,Dsum_lat_reent,~,~] = ariane_Dop(init_volume_g,D);
        Dsum_frac_reent = Dsum_lat_reent/sum(init_volume_g(S)); % Build fraction functionality into ariane_Dop
        Dsum_frac_reent(isnan(Dsum_frac_reent))=0;
	% Density shift
	S = ariane_S({final_age_g,final_lat_g,final_section_g,init_x_g,init_y_g},{[31536000/4 Inf],[-Inf 67],7,[300 Inf],[440 Inf]});
	edges_rho = 27.4:0.05:28.6;
	Dinit = ariane_D(init_dens_g,edges_rho,S);
	[~,Dinitsum,~,~] = ariane_Dop(init_volume_g, Dinit);
	Dfinal = ariane_D(final_dens_g,edges_rho,S);
	[~,Dfinalsum,~,~] = ariane_Dop(init_volume_g, Dfinal);
	
	resetcalcs
	donegincalcs='y'
end
% FIGURE
% Figure properties
gaph = 0.04; gapv = 0.08;
llat = 0.1; blat = 0.5; wlat = 0.2; hlat = 0.35;
lmap = llat+wlat+gaph; bmap = blat+gapv/2; wmap = wlat*2.5; hmap = hlat;
hrho = hmap/2; lrho = llat; brho = blat-gapv-hrho; wrho = wmap+gaph+wlat;

lcbarmap = lmap+gaph; bcbarmap = blat; wcbarmap = wmap-2*gaph; hcbarmap = gapv/4;

posmap = [lmap bmap wmap hmap];
poslat = [llat blat wlat hlat];
posrho = [lrho brho wrho hrho];

poscbarmap = [lcbarmap bcbarmap wcbarmap hcbarmap];

cmrho = cbrewer('qual','Paired',12);

% Figure
close all
f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);

% Latitude
axes('position',poslat);
Dsum_tot = cumsum(Dsum_frac);
hold on;
barh(centre(edges_lat),Dsum_tot+Dsum_frac_reent,1,'facecolor',cmrho(2,:),'edgecolor','none');
barh(centre(edges_lat),Dsum_tot,1,'facecolor',cmrho(1,:));
ylim(latlimB);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
ylabel('Latitude ($^\circ N$)'); xlabel('Fraction');
grid on;
box on;
title('(a) Final latitude','fontsize',titlefs)

% Re-entrainment
axes('position',posmap);
worldmap(latlimN,lonlimN);
pcolorm(double(gphit_g),double(glamt_g),Dsum_map*1e-11);
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none');
colormap(cment); caxis(climsent)
title('(b) Locations of re-entrainment','fontsize',titlefs)
% Colorbar
xc = linspace(climsent(1),climsent(2),ncent);
yc = 0:2;
[xcg,~]=meshgrid(xc,yc);
zc=xcg;

axes('position',poscbarmap)
contourf(xc,yc,zc,xc,'LineColor','none');
colormap(gca,cment);
set(gca,'ytick',[],'xaxislocation','bottom');
xlabel('Volume ($10^{11}m^3$)','fontsize',axesfs)

% Density
axes('position',posrho);
bar(centre(edges_rho),Dinitsum*1e-14,1,'facecolor',cmrho(4,:));
hold on
bar(centre(edges_rho),Dfinalsum*1e-14,1,'facecolor',cmrho(4,:),'facealpha',0.5);
xlim([edges_rho(1), edges_rho(end)]);
xlabel('Neutral density ($kgm^{-3}$)'); ylabel('Volume ($10^{14}m^3$)');
title('(c) Density distribution of water re-entrained outside the GIN seas','fontsize',titlefs)
l = legend('Initial density','Density at re-entrainment','Location','northwest');
set(l,'interpreter','latex','box','off');

if strcmp(savefigs,'y')
%    export_fig([savedir 'fig_GIN'],'-png',saveres);
%	saveas(f,[savedir 'fig_GIN'],'png');
	print(f, [savedir 'fig_GIN'], '-dpng',saveres);
end

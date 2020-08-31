close all
universal
savefigs='n';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBDUCTION LOCATION
% CALCULATIONS
if donesubcalcs=='n';
	if exist('matfiles/D2_fXfY_sum_volume.mat')~=2;
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

	resetcalcs	
	donesubcalcs='y';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCATION OF SUBDUCTION + WINTER-MEAN MLD + WINTER-MEAN PSI

% Figure properties
gapv = 0.1; gaph = 0.05;
lc = 0.1; bc = 0.2; wc = 0.4; hc = 0.2;
lcbarc = lc+wc-gaph; bcbarc = bc+hc/2; wcbarc = 0.01; hcbarc = hc/2;
lb = lc+wc+gaph; bb = bc; wb = wc; hb = hc;
lcbarb = lb+wb-gaph; bcbarb = bb+hb/2; wcbarb = wcbarc; hcbarb = hb/2;
la = lc; ba = bb+hb+gapv; wa = 2*wc; ha = 2*hc;
lcbara = la+wa-2*gaph; bcbara = ba+ha/2; wcbara = 2*wcbarc; hcbara = ha/2;

possub = [la ba wa ha];
posmld = [lb bb wb hb];
pospsi = [lc bc wc hc];
poscbarsub = [lcbara bcbara wcbara hcbara];
poscbarmld = [lcbarb bcbarb wcbarb hcbarb];
poscbarpsi = [lcbarc bcbarc wcbarc hcbarc];

f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 9.7 9.7]);

%%% Location of subduction
axes('position',possub)
worldmap(latlimN,lonlimN);
pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),Dsum*1e-13); shading flat
contourm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),double(mbathy(TrangeX,TrangeY)),bathys,'color',bathycolor)
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmsub); caxis(climssub)

title('(a) Volume subducted','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(climssub(1),climssub(2),ncsub);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarsub)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cmsub);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('$10^{13}\,m^3$','fontsize',axesfs)

%%% Winter-mean mixed layer depth
axes('position',posmld)
worldmap(latlimN,lonlimN);
pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),somxl010_wmean(TrangeX,TrangeY)); shading flat
contourm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),double(mbathy(TrangeX,TrangeY)),bathys,'color',bathycolor)
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmmld); caxis(climsmld)

title('(c) Mixed layer depth (late winter mean)','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(climsmld(1),climsmld(2),ncmld);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarmld)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cmmld);
set(gca,'xtick',[],'ydir','r','yaxislocation','right');
ylabel('$m$','fontsize',axesfs)

%%% Winter-mean Barotropic streamfunction
sobarstf_wmean(sobarstf_wmean>0)=NaN;

axes('position',pospsi)
worldmap(latlimN,lonlimN);
contourm(double(gphit(TrangeX,TrangeY)),...
    double(glamt(TrangeX,TrangeY)),...
    sobarstf_wmean(TrangeX,TrangeY)*1e-6,...
    linspace(climspsi(1),climspsi(2),ncpsi/2));
contourm(double(gphit(TrangeX,TrangeY)),...
    double(glamt(TrangeX,TrangeY)),...
    double(mbathy(TrangeX,TrangeY)),bathys,'color',bathycolor)
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmpsi); caxis(climspsi)

title('(b) Barotropic streamfunction, $\psi$ (late winter mean)','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(climspsi(1),climspsi(2),ncpsi);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarpsi)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cmpsi);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('$Sv$','fontsize',axesfs)

if strcmp(savefigs,'y')
	%export_fig([savedir 'fig_sub-mld-psi'],'-png',saveres);
%	saveas(f,[savedir 'fig_sub-mld-psi'],'png');
	print(f, [savedir 'fig_sub-mld-psi'], '-dpng',saveres);
end

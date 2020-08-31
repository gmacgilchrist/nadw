close all
universal
savefigs='y';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

if donediamapcalcs=='n';
	% Load the 2D histogram of mean density changes
	load(filedia);
	load(filedia_onlyD);
	% Load time-mean kz
	votkeavt = ncread(filew,'votkeavt');
	mbathy  = ncread([rootdir 'grids_ncmod/' gridzf ],'mbathy');
       	[nx,ny,nz] = size(votkeavt);
	kz = nan(nx,ny);
	for i=1:nx;
    		for j=1:ny;
       			if mbathy(i,j)==0;
            			kz(i,j)=NaN;
        		else
            			kz(i,j) = votkeavt(i,j,mbathy(i,j)-1);
        		end
    		end
	end	
	% Load 2D histogram for mean kz along pathway
	load(filekz);
end

gaph = 0.1; gapv = 0.1;
lc = 0.1; bc = 0.1; wc = 0.8; hc = 0.2;
la = lc; ba = bc+hc+gapv; wa = wc/2 - gaph/2; ha = hc*2;
lb = la+wa+gaph; bb = ba; wb = wa; hb = ha;

chtot = gapv/6;
lcaa = la+gaph/4; bcaa = ba-chtot/2; wcaa = wa-gaph/2; hcaa = chtot/2;
lcab = lcaa; bcab = bcaa+hcaa; wcab = wcaa; hcab = hcaa;

lcb = lb+gaph/4; bcb = bcaa; wcb = wcaa; hcb = chtot;
lcc = lc+wc-gaph/4; bcc = bc; wcc = gaph/4; hcc = hc;

posdrdt = [la ba wa ha];
poskzmap = [lb bb wb hb];
poshist = [lc bc wc-gaph/2 hc];

poscbaraa = [lcaa bcaa wcaa hcaa];
poscbarab = [lcab bcab wcab hcab];
poscbarb = [lcb bcb wcb hcb];
poscbarc = [lcc bcc wcc hcc];

f1 = figure('units','inches');
if strcmp(savefigs,'y');
	set(f1,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);

m = 10; r=4;
cm = cbrewer('div','RdBu',2*r*10);
% Blue red color bar
cm_new = [flipud(cm(1:size(cm,1)/2,:));...
    repmat(cm(1,:),[(m-r)*10 1]).*ones((m-r)*10,3);...
    repmat(cm(end,:),[(m-r)*10 1]).*ones((m-r)*10,3);...
    flipud(cm(size(cm,1)/2+1:end,:))];

axes('position',posdrdt)
worldmap(latlimNH,lonlimNH);
% Adjust range of saved data to match with Trange
[nx,ny]=size(Dmean);
rangeX = 4:nx-22;
rangeY = 3:ny-1;
% Mask out areas with less than 10 passages
Dmean(D.n<1000)=nan;
% Plot negative log of positive values
% Log returns negative values, so positive -> positive = blue in figure ->
% getting denser
Dmean_pos = Dmean(rangeX,rangeY);
Dmean_pos(Dmean_pos<0)=NaN;
pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),-log10(Dmean_pos));
shading flat
hold on
% Plot log of absolute negative values
% Log returns negative values, so negative -> negative = red in figure ->
% getting lighter
Dmean_neg = Dmean(rangeX,rangeY);
Dmean_neg(Dmean_neg>0)=NaN;
p=pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),log10(abs(Dmean_neg)));
shading flat
% Set nans to transparent
p.FaceAlpha='texturemap';
alpha(p,double(~isnan(Dmean_neg)))
caxis([-m m])
colormap(cm_new)

contourm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),double(mbathy(TrangeX,TrangeY)),[bathys(4) bathys(4)],'color',bathycolor)
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');

title('(a) Rate of change of density','fontsize',titlefs)

% colorbar
yc = 0:2;
xc = linspace(-m,-m+r,r*10);
[xcg,ycg]=meshgrid(xc,yc);
zc=xcg;

axc1 = axes('position',poscbarab);
contourf(xc,yc,zc,xc,'LineColor','none');
colormap(axc1,flipud(cm(1:size(cm,1)/2,:)));
set(axc1,'ytick',[],'yticklabel',[],'xticklabel',[],'box','on')
axc2 = axes('position',poscbaraa);
contourf(xc,yc,zc,xc,'LineColor','none');
colormap(axc2,cm(size(cm,1)/2+1:end,:));
set(gca,'ytick',[],'xaxislocation','bottom','box','off');
xlabel(axc2,'$\Delta\gamma^n/\Delta t\,[kg\,m^{-3}\,s^{-1}]$','interpreter','latex')

%%%%%
% Map of Kz
axes('position',poskzmap);
worldmap(latlimNH,lonlimNH)
pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),log10(kz(TrangeX,TrangeY)));
shading flat
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');

colormap(gca,cmkz); caxis(climskz)

contourm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),double(mbathy(TrangeX,TrangeY)),[bathys(4) bathys(4)],'color',bathycolor);

title('(b) Coefficient of vertical diffusion','fontsize',titlefs)
% Colorbar
yc = 0:2;
xc = linspace(climskz(1),climskz(2),nckz);
[xcg,ycg]=meshgrid(xc,yc);
zc=xcg;

axes('position',poscbarb)
contourf(xc,yc,zc,xc,'LineColor','none');
colormap(gca,cmkz);
set(gca,'ytick',[],'xaxislocation','bottom');
xlabel('$log_{10}(k_z)$','fontsize',axesfs)

if strcmp(savefigs,'y')
%       export_fig([savedir 'fig_diapycnal_map.png'],'-png',saveres);
%        saveas(f1,[savedir 'fig_diapycnal_map'],'png');
	print(f1, [savedir 'fig_diapycnal_map'], '-dpng',saveres);
end

f2 = figure('units','inches');
if strcmp(savefigs,'y');
        set(f2,'visible','off');
end
%%%%%
% Histogram
gaph = 0.01;
l = 0.1; b = 0.1; h = 0.8; w = 0.3;
lc = l+w-gaph-gaph*2; bc = b+2*gaph; hc = h/4; wc = gaph*2;
poshist = [l b w h];
poscbarc = [lc bc wc hc];
axes('position',poshist);
ntot = nansum(Dsum,1);
Dfrac = Dsum./repmat(ntot,[size(Dsum,1) 1]);
imagesc(centre(D.edges1),centre(D.edges2),Dfrac');
set(gca,'ydir','r');
grid minor
xlim([-5.2 -2]);
clims = [0 0.5]; nc = 50; cm = cmocean('-ice',nc);
colormap(gca,cm); caxis(clims);
ylabel('Neutral density, $\gamma^n\,[kgm^{-3}]$')
xlabel('$log_{10}(k_z)$');
title('Vertical diffusion in NADW','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(clims(1),clims(2),nc);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarc)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cm);
set(gca,'xtick',[],'yaxislocation','left');
ylabel('Fraction','fontsize',axesfs)


if strcmp(savefigs,'y')

%   	export_fig([savedir 'fig_diapycnal_map.png'],'-png',saveres);
%	saveas(f2,[savedir 'fig_diapycnal_kz'],'png');
	print(f2, [savedir 'fig_diapycnal_kz'], '-dpng',saveres);
end

% DRAW FIGURES FOR JOURNAL OF CLIMATE PAPER
%
clear; close all;

savefigs = 'n';
savedir = '/home/ocean2/graemem/frontend/paperfigures/2018_nadw/';
saveres = '-r90';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
% file locations
rootdir = '/home/ocean2/graemem/ariane/';
model = 'orca025_global_5d';
experiment = 'quant_back_seedNAn1_t3560-sep-4217_sign27.7-28_MLrefz8delsig0.01';
fileout = [rootdir 'experiments/' model '/' experiment '/ariane_positions_quantitative.nc'];
filein = [rootdir 'experiments/' model '/' experiment '/ariane_initial.nc'];
gridf = 'ORCA025.L75-GRD88_mesh_hgr.nc';

% model data
load([rootdir 'time/time_' model '.mat'],'time');
tmask = ncread([rootdir 'grids_ncmod/' gridf ],'tmask');
gphit = ncread([rootdir 'grids_ncmod/' gridf ],'gphit');
glamt = ncread([rootdir 'grids_ncmod/' gridf ],'glamt');
% model data direct from model
modeldir = '/home/ocean_shared_data1/DRAKKAR/ORCA025.L75-GJM189-S/';
meanfile = 'MarchApril_means_1975-2015/ORCA025.L75-GJM189_Mar-Apr-1975-2015-mean_';
somxl010_wmean = ncread([modeldir meanfile 'gridT.nc'],'somxl010');
sobarstf_wmean = ncread([modeldir meanfile 'psi.nc'],'sobarstf');
gridfile = 'GRID/ORCA025.L75-GJM189_mesh_zgr.nc';
mbathy = ncread([modeldir gridfile],'mbathy');
gdept = ncread([modeldir gridfile],'gdept_0');

% ariane data
final_section = ncread(fileout,'final_section');
init_t = ncread(fileout,'init_t');
final_age = ncread(fileout,'final_age');
init_temp = ncread(fileout,'init_temp');
final_temp = ncread(fileout,'final_temp');
init_salt = ncread(fileout,'init_salt');
final_salt = ncread(fileout,'final_salt');
init_dens = ncread(fileout,'init_dens');
init_volume = ncread(filein,'init_volume');
init_x = ncread(fileout,'init_x');
init_y = ncread(fileout,'init_y');
init_z = ncread(fileout,'init_z');
final_x = ncread(fileout,'final_x');
final_y = ncread(fileout,'final_y');
init_lat = ncread(fileout,'init_lat');
final_lat = ncread(fileout,'final_lat');

region_limits = load([rootdir 'experiments/' model '/' experiment '/region_limits']);

% Neutral density at time of subduction
final_dens = calc_sigmantr(final_temp,final_salt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPECIFICATIONS FOR ALL FIGURES

% Mask colors
landcolor = 0.6*[1 1 1];
unventcolor = 0.4*[1 1 1];

% Basin-wide map limits
latlimB = [-30 80]; lonlimB = [-100 20];
% Northern region map limits
latlimN = [45 75]; lonlimN = [-65 10];
% Labrador Sea map limits
latlimLS = [48 65]; lonlimLS = [-62 -42];

% Bathymetric contours to show
bathys = [40 51 60 65];
bathycolor=0.8*[1 1 1];

% Font sizes
titlefs = 12;
axesfs = 10;

% Colormaps
% Subduction volume
climssub = [0 3]; ncsub = 50; cmsub = cmocean('dense',ncsub);
% Mixed layer depth
climsmld = [0 1500]; ncmld = 50; cmmld = cmocean('deep',ncmld);
% Barotropic streamfunction
climspsi = [-40 -0]; ncpsi = 40; cmpsi = cmocean('thermal',ncpsi);
% Re-entrainment volume
climsent = [0 5]; ncent = 50; cment = cmocean('dense',ncent);
% Age
climsage = [0 58]; ncage = 58; cmage = cmocean('matter',ncage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL CALCULATIONS

% Determine YEAR of ventilation
[final_year,final_agereal] = calc_yr_ariane(final_age,init_t,model,time);
% Determine MONTH of ventilation
mtdec = (final_year-floor(final_year)); % Months as a decimal fraction of 1
        
% Put into discrete months
fm = 4; % fraction of month to divide into
final_month = nan(size(mtdec));
for mt = 1:12*fm;
    final_month(mtdec>=(mt-1)/(12*fm) & mtdec < mt/(12*fm))=mt/fm;
end;

% Find the unseeded grid squares
S = ariane_S({init_t},{4217});
ind = sub2ind(size(tmask),ceil(init_x),ceil(init_y),floor(init_z));
section = nan(size(tmask));
section(ind(S)) = final_section(S);
seeded = sum(isfinite(section),3);

% Find ventilated particles at last timestep for qualitative experiment
ind = find(init_t == 4217 & final_section==7 & init_dens>27.7 & init_dens<27.8);
fileid = fopen([rootdir 'experiments/' model '/qual' experiment(6:end) '_subbin-t4217-sec7-sig27.7-27.8/subset.txt'],'w');
fprintf(fileid,'%i \n',ind);
fclose(fileid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MONTH OF VENTILATION
% CALCULATIONS
% Lagrangian
S = ariane_S({init_t,final_section,final_age},{4217,7,[31536000 Inf]});
edges = 0+1/(fm*2):1/fm:12+1/(fm*2);
D = ariane_D(final_month,edges,S);
[~,Dsum,~,~] = ariane_Dop(init_volume,D);
% Eulerian
years = 1958:2015;
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
edges_eul = edges/12;
D = histcounts(fractNA,edges_eul);
%%
close all
% Figure properties
lmap = 0.1; bmap = 0.5; wmap = 0.5; hmap = 0.4;
gapv = 0.03;
hhist = 0.1; lhist = lmap;  bhist = bmap-gapv-hhist; whist = wmap;
hcbar = 0.01; lcbar = lmap; bcbar = bhist-gapv/2-hcbar; wcbar = wmap;
posmap = [lmap bmap wmap hmap];
poshist = [lhist bhist whist hhist];
poscbar = [lcbar bcbar wcbar hcbar];
clims = [0 12]; nc = 12; cm = cbrewer('qual','Paired',nc);
cm = [unventcolor; cm];
    
figure('units','inches');
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

title('(a) Month of subduction in NADW (depth-median)','fontsize',titlefs)

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
ylabel('Volume ($10^{15}\, m^3$)');
% Eulerian
ax_eul = axes('position',poshist);
for mt = 1:12;
    fc = cmc(mt,:);
    if mt==12
        ind = (mt-1)*fm+1:(mt-1)*fm+fm;
    else
        ind = (mt-1)*fm+1:(mt-1)*fm+fm+1;
    end
    plot(edges_eul(ind),D(ind)/sum(D),'color','k','linewidth',2)
    hold on
end
set(ax_eul,'color','none','yaxislocation','right','xticklabels',[])
ylabel('Fraction')
title('(b) Volume of NADW subducted in each month','fontsize',titlefs)

% Colour bar
yc = 0.5:12.5;
xc = 0:2;
[~,zc] = meshgrid(xc,yc);
ax = axes('position',poscbar);
contourf(yc,xc,zc',yc,'linestyle','none'); caxis(clims); colormap(ax,cbrewer('qual','Paired',12));
set(gca,'ytick',[],'xtick',1:12,'ticklength',[0 0],...
    'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
xlabel('Month of subduction')

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_month'],'-png',saveres);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBDUCTION LOCATION
% CALCULATIONS
S = ariane_S({init_t,final_section,final_age},{4217,7,[31536000 Inf]});

% Distribution at subduction location
minX = floor(min(final_x(S))); maxX = ceil(max(final_x(S)));
minY = floor(min(final_y(S))); maxY = ceil(max(final_y(S)));
edgesX = minX:1:maxX;
edgesY = minY:1:maxY;
D = ariane_D2(final_x,final_y,edgesX,edgesY,S);
[Dind,Dsum,~,~] = ariane_Dop(init_volume,D);
Dsum(Dsum==0)=NaN;
% Range of tracer points on model grid (shifted right and up for
% correspondence with model indexing)
TrangeX = minX+1:maxX; TrangeY = minY+1:maxY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCATION OF SUBDUCTION + WINTER-MEAN MLD + WINTER-MEAN PSI
close all

% Figure properties
gapv = 0.075; gaph = 0;
lc = 0.1; bc = 0.05; wc = 0.6; hc = 0.25;
lcbarc = lc+wc-gaph; bcbarc = bc+hc/2; wcbarc = 0.02; hcbarc = hc/2;
lb = lc+wc+gaph; bb = bc; wb = wc; hb = hc;
lcbarb = lb+wb-gaph; bcbarb = bb+hb/2; wcbarb = wcbarc; hcbarb = hb/2;
la = lc; ba = bb+hb+gapv; wa = 2*wc; ha = 2*hc;
lcbara = la+wa-gaph; bcbara = ba+ha/2; wcbara = wcbarc; hcbara = ha/2;

possub = [la ba wa ha];
posmld = [lb bb wb hb];
pospsi = [lc bc wc hc];
poscbarsub = [lcbara bcbara wcbara hcbara];
poscbarmld = [lcbarb bcbarb wcbarb hcbarb];
poscbarpsi = [lcbarc bcbarc wcbarc hcbarc];

figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);

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

title('(b) Mixed layer depth (late winter mean)','fontsize',titlefs)

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
contourfm(double(gphit(TrangeX,TrangeY)),...
    double(glamt(TrangeX,TrangeY)),...
    sobarstf_wmean(TrangeX,TrangeY)*1e-6,...
    linspace(climspsi(1),climspsi(2),ncpsi/2),'edgecolor','none');
contourm(double(gphit(TrangeX,TrangeY)),...
    double(glamt(TrangeX,TrangeY)),...
    double(mbathy(TrangeX,TrangeY)),bathys,'color',bathycolor)
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmpsi); caxis(climspsi)

title('(c) Barotropic streamfunction, $\psi$ (late winter mean)','fontsize',titlefs)

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
    export_fig([savedir 'fig_sub-mld-psi'],'-png',saveres);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VENTILATION AS FUNCTION OF BATHYMETRY AND DENSITY
%
close all
% BATHYMETRIC CONTOUR OF SUBDUCTION
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

figure
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
title('(a) Bathymetric model level over which subduction occurs','fontsize',titlefs)

% DENSITY DISTRIBUTION IN EACH BATHYMETRY BAND
cm = cbrewer('qual','Set1',9);

rr = 27.7:0.01:28;
% Initial density
Dbid = ariane_D2(final_bathy,init_dens(S),[-Inf bathys Inf],rr,Sb);
[~,Dbidsum,~,~] = ariane_Dop(init_volume(S),Dbid);
Dbfd = ariane_D2(final_bathy,final_dens(S),[-Inf bathys Inf],rr,Sb);
[~,Dbfdsum,~,~] = ariane_Dop(init_volume(S),Dbfd);
subplot(212); hold on
h=cell(3,1);
s=cell(3,1);
for b=2:4;
    pid = Dbidsum(b,:);
    h{b-1}=plot(centre(rr),pid*1e-14,'-','color',cm(b,:),'linewidth',2);
    pfd = Dbfdsum(b,:);
    plot(centre(rr),pfd*1e-14,'--','color',cm(b,:));
    s{b-1}=[num2str(round(gdept(bathys(b-1)),-1)) ' to ' num2str(round(gdept(bathys(b)),-1)) ' m'];
end
ylims = [0 8];
set(gca,'box','on','ylim',ylims)
xlabel('Neutral density, $\gamma^n$'); ylabel('Volume ($10^{14}m^3$)');
title('(b) Density distribution for subduction over bathymetric ranges','fontsize',titlefs)
l = legend([h{1} h{2} h{3}],{s{1},s{2},s{3}});
set(l,'box','off','location','northeast','interpreter','latex')

plot([27.71 27.73],[ylims(2)-1 ylims(2)-1],'k--'); text(27.735,ylims(2)-1,'Density at subduction');
plot([27.71 27.73],[ylims(2)-1.8 ylims(2)-1.8],'k-','linewidth',2); text(27.735,ylims(2)-1.8,'Density in Sep. 2015');

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_subduction_bathy-dens'],'-png',saveres);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MLD ALONG PSI, SUBDUCTION VOLUME ALONG PSI
close all;
% FINDING AND PLOTTING CONTOURS
psir = -40:1:0; % Find contours in certain range
[C,h] = contour(sobarstf_wmean(TrangeX,TrangeY)*1e-6,psir);
close
start = 1; stop = start+1+C(2,start)-1; l = 1;
cpos = cell(1); clev = nan(1);
while stop<length(C);
    stop = start+1+C(2,start)-1;
    cpos{l} = C(:,start+1:stop);
    clev(l) = C(1,start);
    start = stop+1;
    l = l+1;
end

close all

% Contours in subpolar gyre
cpos_subp = cell(1);
clev_subp = cell(1);
n = 1;
for l = psir;
    ind = find(clev==l);
    for c=ind;
        % Determine if the contours are present in a region of the subpolar gyre
        yes=0;
        for i=1:length(cpos{c}(1,:));
            if cpos{c}(1,i)>410 && cpos{c}(1,i)<450 && cpos{c}(2,i)>160 && cpos{c}(2,i)<220;
                yes=1;
            end
        end
        if yes==1;
            cpos_subp{n} = cpos{c};
            clev_subp{n} = l;
            n = n+1;
        end
    end
end

% Mixed layer depth along contours
% truncated fields to correspond to subducted volume distribution
mld = somxl010_wmean(TrangeX,TrangeY);
lat = gphit(TrangeX,TrangeY);
lon = glamt(TrangeX,TrangeY);
% preallocate cell arrays
cmld = cell(size(cpos_subp));
clat = cell(size(cpos_subp));
clon = cell(size(cpos_subp));
cDsum = cell(size(cpos_subp));
for c = 1:length(cpos_subp);
    % x and y positions (note index)
    xs = round(cpos_subp{c}(2,:));
    ys = round(cpos_subp{c}(1,:));
    
    for i=1:length(xs);
        cmld{c}(i) = mld(xs(i),ys(i));
        clat{c}(i) = lat(xs(i),ys(i));
        clon{c}(i) = lon(xs(i),ys(i));
        cDsum{c}(i) = Dsum(xs(i),ys(i));
    end
end

% Plot every contour from 1 to 53, and identify ones in boundary current
% and central Labrador Sea
% close all
% figure;
% worldmap([50 65],[-65 -40]);
% geoshow('landareas.shp'); hold on;
% for c=1:length(cmld)
%     plotm(double(clat{c}),double(clon{c}),'-');
%     title(num2str(c));
%     pause
% end

conts_bc = [30 32 34 36:53]; % confined to bathymetry
conts_rcopen = [18 19 22 24 26 27]; % recirculate within LS, also extend into subP gyre
conts_rccloseW = [1:5 7 10 15]; % closed contours in W recirc
conts_rccloseE = [31 33 35]; % closed contours in E recirc

% close all
% figure;
% worldmap([50 65],[-65 -40]);
% geoshow('landareas.shp'); hold on;
% for c=conts_bc;
%     plotm(double(clat{c}),double(clon{c}),'r-');
% end
% for c=conts_rcopen;
%     plotm(double(clat{c}),double(clon{c}),'b-');
% end
% for c=conts_rccloseW;
%     plotm(double(clat{c}),double(clon{c}),'g-');
% end
% for c=conts_rccloseE;
%     plotm(double(clat{c}),double(clon{c}),'k-');
% end

%% FINAL PLOTS
close all
% PLOT MLD ALONG THE BC COUNTOURS
% Figure properties
gaph = 0.05; gapv = 0.01;
lmap = 0.1; bmap = 0.5; wmap = 0.4; hmap = 0.4;
lsub = lmap+wmap+gaph; bsub = bmap; wsub = wmap; hsub = hmap/2-2*gapv;
lmld = lsub; bmld = bsub+hsub+gapv; wmld = wsub; hmld = hsub;

hcbarmld = 0.01;
lcbarmld = lmap; bcbarmld = bmap-4*hcbarmld-gapv; wcbarmld = wmap;
lcbarpsi = lsub; bcbarpsi = bcbarmld; wcbarpsi = wsub; hcbarpsi = hcbarmld;

posmap = [lmap bmap wmap hmap];
possub = [lsub bsub wsub hsub];
posmld = [lmld bmld wmld hmld];

poscbarmld = [lcbarmld bcbarmld wcbarmld hcbarmld];
poscbarpsi = [lcbarpsi bcbarpsi wcbarpsi hcbarpsi];

% define a starting box
stalon=[-45 -43]; stalat=[58 60];
stolon=[-47 -46]; stolat=[46 51];

figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);
axes('position',posmap);
axmap = worldmap(latlimLS,lonlimLS);
axmld = axes('position',posmld); hold(axmld,'on')
axsub = axes('position',possub); hold(axsub,'on')
for c = conts_bc;
    sta=find(clon{c}>stalon(1) & clon{c}<stalon(2) & clat{c}>stalat(1) & clat{c}<stalat(2),1,'first');
    sto=find(clon{c}>stolon(1) & clon{c}<stolon(2) & clat{c}>stolat(1) & clat{c}<stolat(2),1,'first');
    if c>=51;
        ind = sta:sto;
    else
        ind = [sta:length(cmld{c}) 1:sto];
    end
    scatterm(axmap,double(clat{c}(ind)),double(clon{c}(ind)),cDsum{c}(ind)*1e-12, cmld{c}(ind),'filled');
    dist = [0 cumsum(sw_dist(clat{c}(ind),clon{c}(ind)))];
    plot(axmld,dist,cmld{c}(ind),'color',cmpsi(40+clev_subp{c},:)); 
    plot(axsub,dist,cDsum{c}(ind)*1e-13,'color',cmpsi(40+clev_subp{c},:));
end
axes(axmap)
textm(57.5,-52.5+0.4,'$10^{13}\,m^3$'); textm(63,-47.5,'(a)','fontsize',titlefs)
scatterm(axmap,57.5,-52.5,10,cmmld(round(ncmld/2),:),'filled');
contourm(double(gphit(TrangeX,TrangeY)),...
    double(glamt(TrangeX,TrangeY)),...
    double(mbathy(TrangeX,TrangeY)),bathys,'color',bathycolor)
geoshow(axmap,'landareas.shp','facecolor',landcolor,'edgecolor','none'); hold on;
colormap(axmap,cmmld); caxis(axmap,climsmld);
set(axmld,'xlim',[0 1700],'box','on','xgrid','on','ygrid','on','xticklabel',[])
ylabel(axmld,'MLD ($m$)'); text(axmld,100,1000,'(b)','fontsize',titlefs,...
    'horizontalalignment','c','verticalalignment','bottom')
set(axsub,'xlim',[0 1700],'box','on','xgrid','on','ygrid','on')
xlabel(axsub,'Distance along contour ($km$)'); ylabel(axsub,'Volume ($10^{13}m^3$)');
text(axsub,100,4,'(c)','fontsize',titlefs,...
    'horizontalalignment','c','verticalalignment','bottom')

% MLD Colorbar
xc = linspace(climsmld(1),climsmld(2),ncmld);
yc = 0:2;
[xcg,~]=meshgrid(xc,yc);
zc=xcg;

axes('position',poscbarmld)
contourf(xc,yc,zc,xc,'LineColor','none');
colormap(gca,cmmld);
set(gca,'ytick',[],'xaxislocation','bottom');
xlabel('Mixed layer depth ($m$)','fontsize',axesfs)

% psi Colorbar
xc = linspace(climspsi(1),climspsi(2),ncpsi);
yc = 0:2;
[xcg,~]=meshgrid(xc,yc);
zc=xcg;

axes('position',poscbarpsi)
contourf(xc,yc,zc,xc,'LineColor','none');
colormap(gca,cmpsi);
set(gca,'ytick',[],'xaxislocation','bottom');
xlabel('$\psi$ ($Sv$)','fontsize',axesfs)

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_mld-on-psi'],'-png',saveres);
end

%% EXPORT PATHWAYS AS A FUNCTION OF BATHYMETRY OF SUBDUCTION
% Retrieve particles subducted in different bathymetric ranges
close all
for b=1:length(bathys)-1;
    % Index of particles in bathy range, relative to subset S
cells = Dbind(bathys(b):bathys(b+1));
indb = [];
for i=1:length(cells);
    indb=[indb cells{i}];
end
disp(length(indb))
Sb = S(S);

% Extract particles in subset S that subducted in bathymetry range
x = init_x(S); y = init_y(S); z = init_z(S);
yr = final_year(S);
indL = sub2ind(size(tmask),ceil(x(indb)),ceil(y(indb)),floor(z(indb))); 
year = nan(size(tmask));
year(indL) = yr(indb);

year_med = nanmedian(year,3);

figure;
worldmap(latlimB,lonlimB);
pcolorm(double(gphit),double(glamt),year_med);
geoshow('landareas.shp')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VENTILATION VARIABILITY
% Load the remaining time-slices
final_age_all = final_age;
final_section_all = final_section;
final_x_all = final_x;
final_y_all = final_y;
init_x_all = init_x;
init_y_all = init_y;
init_z_all = init_z;
final_dens_all = final_dens;
init_dens_all = init_dens;
init_t_all = init_t;
init_lat_all = init_lat;
final_lat_all = final_lat;
init_volume_all = init_volume;
timeslices = {'2829-sep-3487','2099-sep-2756','1368-sep-2026'};
for t = 1:length(timeslices);
    experiment = ['quant_back_seedNAn1_t' timeslices{t} '_sign27.7-28_MLrefz8delsig0.01'];
    fileout = [rootdir 'experiments/' model '/' experiment '/ariane_positions_quantitative.nc'];
    filein = [rootdir 'experiments/' model '/' experiment '/ariane_initial.nc'];
    % Tack them on
    final_age_all = [ncread(fileout,'final_age'); final_age_all];
    final_x_all = [ncread(fileout,'final_x'); final_x_all];
    final_y_all = [ncread(fileout,'final_y'); final_y_all];
    init_x_all = [ncread(fileout,'init_x'); init_x_all];
    init_y_all = [ncread(fileout,'init_y'); init_y_all];
    init_z_all = [ncread(fileout,'init_z'); init_z_all];
    final_section_all = [ncread(fileout,'final_section'); final_section_all];
    final_dens_all = [calc_sigmantr(ncread(fileout,'final_temp'),ncread(fileout,'final_salt')); final_dens_all];
    init_dens_all = [ncread(fileout,'init_dens'); init_dens_all];
    init_t_all = [ncread(fileout,'init_t'); init_t_all];
    init_lat_all = [ncread(fileout,'init_lat'); init_lat_all];
    final_lat_all = [ncread(fileout,'final_lat'); final_lat_all];
    init_volume_all = [ncread(filein,'init_volume'); init_volume_all];
end
times = unique(init_t_all);
% Load the calculated bathymetry-age distribution
load Dsum_agebathy
load Dsum_agedens
load Dsum_age
% Load water mass transformation
load /home/ocean_personal_data/graemem/dwfproxy/DRAKKAR/matfiles/WMT_deepocean_5dglobal_k8_ds0.01.mat

%%% CALCULATIONS
% Derive distributions one timeslice at a time
edges_age = -3/12:1:57+9/12;
edges_bathy = [-Inf bathys Inf];
edges_dens = 27.7:0.01:28;
years = 1958:2016;
ages = years-1958;
years_seed = 1976:2015;

% AGE ANOMALY
edges_t = [-Inf; centre(times); Inf];
S = ariane_S({final_section_all},{7});
D = ariane_D2(init_t_all,final_age_all/31536000,edges_t,edges_age,S);
[~,Dsum,~,~] = ariane_Dop(init_volume_all,D);

Dnum = sum(isfinite(Dsum_age));
Dsummean_age = nansum(Dsum_age)./Dnum;
Danom_age = (Dsum_age(end,:)-Dsummean_age)./Dsummean_age;

% AGE-BATHYMETRY
for t = 1:length(times);
disp(['Timeslice: ' num2str(times(t)) ]);
S = ariane_S({init_t_all,final_section_all},{times(t),7});

% Locations of subduction
minX = floor(min(final_x_all(S))); maxX = ceil(max(final_x_all(S)));
minY = floor(min(final_y_all(S))); maxY = ceil(max(final_y_all(S)));
edgesX = minX:1:maxX;
edgesY = minY:1:maxY;
D_xy = ariane_D2(final_x_all,final_y_all,edgesX,edgesY,S);
[Dind,Dsum,~,~] = ariane_Dop(init_volume_all,D_xy);
% Range of tracer points on model grid (shifted right and up for
% correspondence with model indexing)
TrangeX = minX+1:maxX; TrangeY = minY+1:maxY;

% Recover bathymetry of subduction (for subset S)
final_bathy = nan(sum(S),1);
[nx,ny]=size(Dind);
for i=1:nx;
    for j=1:ny;
        if ~isempty(Dind{i,j})
            final_bathy(Dind{i,j}(:))=mbathy(TrangeX(i),TrangeY(j));
        end
    end
end
Sb = S(S);

D_agebathy = ariane_D2(final_age_all(S)/31536000,final_bathy,edges_age,edges_bathy,Sb);
[~,Dsum_agebathy(:,:,t),~,~] = ariane_Dop(init_volume_all(S),D_agebathy);
end

% AGE-DENSITY
for t = 1:length(times);
disp(['Timeslice: ' num2str(times(t)) ]);
S = ariane_S({init_t_all,final_section_all},{times(t),7});

D_agedensf = ariane_D2(final_age_all/31536000,final_dens_all,edges_age,edges_dens,S);
[~,Dsum_agedensf(:,:,t),~,~] = ariane_Dop(init_volume_all,D_agedensf);
end

%%% FIGURES
close all
% Figure properties
cm = cbrewer('qual','Set1',length(edges_bathy)-1);

gapv = 0.1;
l1 = 0.1; b1 = 0.4; w1 = 0.8; h1 = 0.2;
l2 = l1; b2 = b1+h1+gapv; w2 = w1; h2 = h1;

% AVERAGE AGE DISTRIBUTION AND ANOMALY
close all
figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);
axes('position',[l1 b1+1*(h1+gapv) w1 h1]);
bar(centre(ages),Dsummean_age*1e-14,1,'facecolor',landcolor);
hold on
plot(centre(ages),Dsum_age(end,:)*1e-14,'k.');
text(56.5,18,'$\bullet$ Age distribution in 2015')
text(0.5,21,'$\bullet$','horizontalalignment','c');
text(0.5,22,num2str(round(Dsum_age(end,1)*1e-14)),'horizontalalignment','c');
set(gca,'ylim',[0 20],'xlim',[0 58],'xdir','r');
grid on
xlabel('Age (years)'); ylabel('Volume ($10^{14}m^3$)');
title('(a) Average age distribution of NADW','fontsize',titlefs);

axes('position',[l1 b1+0*(h1+gapv) w1 h1]);
bar(fliplr(centre(years)),Danom_age,1,'facecolor',landcolor);
set(gca,'ylim',[-1.1 1.1],'xlim',[1958 2016]);
grid on
xlabel('Year'); ylabel('Normalised anomaly');
title('(b) Anomaly in age distribution of NADW in 2015','fontsize',titlefs);

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_ageanom'],'-png',saveres);
end

% SUBDUCTION AND RE-ENTRAINMENT ANOMALIES
close all;
figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);
d1 = Dsum_age(:,1)-nanmean(Dsum_age(:,1));
d2 = Dsum_age(:,2)-nanmean(Dsum_age(:,2));
axes('position',[l1 b1+1*(h1+gapv) w1 h1]);
bar(years_seed+0.5,d1*1e-14,1,'facecolor',0*[1 1 1]);
hold on;
bar(years_seed+0.5-1,d2*1e-14,1,'facecolor',1*[1 1 1],'facealpha',0.9);
set(gca,'xlim',[1958 2016]);
grid on;
xlabel('Year'); ylabel('Volume anomaly ($10^{14}m^3$)');
title('(a) Volume anomaly following subduction and after re-entrainment','fontsize',titlefs);
l = legend('Following subduction','After re-entrainment','location','northwest');
set(l,'interpreter','latex','box','off')

% Effectiveness of Stommel's demon
d1 = 100*(Dsum_age(1:end-1,1)-Dsum_age(2:end,2))./Dsum_age(1:end-1,1); % First-year demon
d2 = 100*(Dsum_age(2:end-1,2)-Dsum_age(3:end,3))./Dsum_age(2:end-1,2); % Second-year demon
axes('position',[l1 b1+0*(h1+gapv) w1 h1]);
bar(years_seed(1:end-1)+0.5,d1,1,'facecolor',0.6*[1 1 1]); hold on;
bar(years_seed(1:end-2)+0.5,d2,1,'facecolor',0.4*[1 1 1]);
plot([years(1) years(end)],[nanmean(d1) nanmean(d1)],'k--');
plot([years(1) years(end)],[nanmean(d2) nanmean(d2)],'k--');
grid on
set(gca,'xlim',[1958 2016],'ydir','r','ylim',[0 100]);
xlabel('Year'); ylabel('\%');
title('(b) Percentage volume reduction after one year and between one and two years','fontsize',titlefs)
l = legend('First-year re-entrainment','Second-year re-entrainment','location','southwest');
set(l,'interpreter','latex','box','off')

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_reent'],'-png',saveres);
end

% AGE DISTRIBUTION BY BATHYMETRY
% Renew figure properties
gapv = 0.08;
l1 = 0.1; b1 = 0.1; w1 = 0.8; h1 = 0.2;

Dsummean_agebathy = nansum(Dsum_agebathy,3)./repmat(Dnum',[1 length(edges_bathy)-1]);
Danom_agebathy = (Dsum_agebathy(:,:,end)-Dsummean_agebathy)./repmat(Dsummean_age',[1 length(edges_bathy)-1]);
Dpos = Danom_agebathy;
Dneg = Danom_agebathy;
Dpos(Dpos<0)=0;
Dneg(Dneg>0)=0;

close all;
figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);
axes('position',[l1 b1+h1+gapv+h1+gapv/4+gapv w1 h1]);
bar(fliplr(centre(years)),Dpos,1,'stacked','edgecolor','none');
l = legend('a','b','c','d','e','location','best');
set(l,'interpreter','latex','box','off')
hold on;
bar(fliplr(centre(years)),Dneg,1,'stacked','edgecolor','none');
plot(fliplr(centre(years)),Danom_age,'k.');
grid on;
set(gca,'xlim',[years(1) years(end)],'ylim',[-1.1 1.1]);
colormap(cm);
xlabel('Year'); ylabel('Normalised anomaly');
title('(a) Anomaly in age distribution, by bathymetry of subduction','fontsize',titlefs)

axes('position',[l1 b1+1*(h1+gapv)+h1/2+gapv/4 w1 h1/2]);
d1 = Dsum_agebathy(1,:,:)-repmat(nanmean(Dsum_agebathy(1,:,:),3),[1 1 Dnum(1)]);
d2 = Dsum_agebathy(2,:,:)-repmat(nanmean(Dsum_agebathy(2,:,:),3),[1 1 Dnum(1)]);
bar(years_seed+0.5,squeeze(d1(:,3,:))*1e-14,1,'facecolor',cm(3,:))
hold on;
bar(years_seed+0.5-1,squeeze(d2(:,3,:))*1e-14,1,'facecolor',cm(3,:)+(1-cm(3,:))*3/4,'facealpha',0.9)
grid on
set(gca,'xlim',[years(1) years(end)],'xticklabel',[],'ylim',[-5 5])
title('(b) Volume anomaly following subduction and after re-entrainment','fontsize',titlefs)
text(1960,2.5,'Boundary current','fontsize',titlefs);
l = legend('Following subduction','After re-entrainment','location','southwest');
set(l,'interpreter','latex','box','off');

axes('position',[l1 b1+1*(h1+gapv) w1 h1/2]);
d1 = Dsum_agebathy(1,:,:)-repmat(nanmean(Dsum_agebathy(1,:,:),3),[1 1 Dnum(1)]);
d2 = Dsum_agebathy(2,:,:)-repmat(nanmean(Dsum_agebathy(2,:,:),3),[1 1 Dnum(1)]);
bar(years_seed+0.5,squeeze(d1(:,4,:))*1e-14,1,'facecolor',cm(4,:))
hold on;
bar(years_seed+0.5-1,squeeze(d2(:,4,:))*1e-14,1,'facecolor',cm(4,:)+(1-cm(4,:))*3/4,'facealpha',0.9)
grid on
set(gca,'xlim',[years(1) years(end)],'ylim',[-10 10])
xlabel('Year'); ylabel('$\quad \quad$ Volume anomaly ($10^{14}m^3$)');
text(1960.5,5,'Open ocean','fontsize',titlefs);
l = legend('Following subduction','After re-entrainment','location','southwest');
set(l,'interpreter','latex','box','off');

% First year demon for bathymetric range
axes('position',[l1 b1+0*(h1+gapv) w1 h1]);
db = 100*(Dsum_agebathy(1,:,1:end-1)-Dsum_agebathy(2,:,2:end))./Dsum_agebathy(1,:,1:end-1);
bar(years_seed(1:end-1)+0.5,squeeze(db(:,4,:)),1,'facecolor',cm(4,:)); hold on;
bar(years_seed(1:end-1)+0.5,squeeze(db(:,3,:)),1,'facecolor',cm(3,:));
plot([years(1) years(end)],[nanmean(db(:,4,:)) nanmean(db(:,4,:))],'--','color',cm(4,:));
plot([years(1) years(end)],[nanmean(db(:,3,:)) nanmean(db(:,3,:))],'--','color',cm(3,:));
grid on;
set(gca,'ydir','r','xlim',[years(1) years(end)])
xlabel('Year'); ylabel('\%');
title('(c) Percentage volume reduction after re-entrainment','fontsize',titlefs)
l = legend('Boundary current','Open ocean','location','northwest');
set(l,'interpreter','latex','box','off');

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_anom-bathy'],'-png',saveres);
end

% AGE DISTRIBUTION BY DENSITY
% Figure properties
gapv = 0.1; gaph = 0.01;
l1 = 0.1; b1 = 0.4; w1 = 0.75; h1 = 0.2;
wc = 0.02;
cm = flipud(cbrewer('div','RdBu',33));

% Calculate Age-density anomalies
Dsummean_agedensi = nansum(Dsum_agedensi,3)./repmat(Dnum',[1 length(edges_dens)-1]);
Danom_agedensi = (Dsum_agedensi(:,:,end)-Dsummean_agedensi)./repmat(Dsummean_age',[1 length(edges_dens)-1]);
Dsummean_agedensf = nansum(Dsum_agedensf,3)./repmat(Dnum',[1 length(edges_dens)-1]);
Danom_agedensf = (Dsum_agedensf(:,:,end)-Dsummean_agedensf)./repmat(Dsummean_age',[1 length(edges_dens)-1]);

% Calculate annual average WMT
Fn = FgW;
Fs = F;
sepind = 1:73:size(Fn,1);
Fn(isnan(Fn))=0;
Fs(isnan(Fs))=0;
Fnm=nan(length(sepind)-1,size(Fn,2));
Fsm=nan(length(sepind)-1,size(Fs,2));
for t=1:length(sepind)-1;
Fnm(t,:)=mean(Fn(sepind(t):sepind(t+1)-1,:),1);
Fsm(t,:)=mean(Fs(sepind(t):sepind(t+1)-1,:),1);
end
Fnmanom = Fnm-repmat(mean(Fnm,1),[size(Fnm,1) 1]);
Fsmanom = Fsm-repmat(mean(Fsm,1),[size(Fsm,1) 1]);

close all
figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);

axes('position',[l1 b1+1*(h1+gapv) w1 h1]);
imagesc(centre(years),edges_dens,fliplr(Danom_agedensi'));
caxis(0.1*[-1 1]);
grid on
xlabel('Year'); ylabel('$\gamma^n$');
title('(a) Anomaly in age distribution, by density','fontsize',titlefs)
t = colorbar;
set(t,'position',[l1+w1+gaph b1+1*(h1+gapv) wc h1],'interpreter','latex')
t.Label.String = 'Normalised anomaly';
t.Label.Interpreter = 'latex';
t.Label.FontSize = axesfs;

axes('position',[l1 b1+0*(h1+gapv) w1 h1]);
imagesc(centre(years),g,Fnmanom'*1e-6);
caxis(30*[-1 1]);
grid on
xlabel('Year'); ylabel('$\gamma^n$');
title('(b) Anomaly in annual-mean, surface-forced water mass transformation','fontsize',titlefs)
t = colorbar;
set(t,'position',[l1+w1+gaph b1+0*(h1+gapv) wc h1]);
t.Label.String = 'Sv';
t.Label.Interpreter = 'latex';
t.Label.FontSize = axesfs;

colormap(cm);

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_anom-dens'],'-png',saveres);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PATHWAYS OF VENTILATION
% AGE OF NADW

age_all = [];
section_all = [];
volume_all = [];
for t = 1:length(times);
    disp(times(t))
    S = ariane_S({init_t_all},{times(t)});
    
    % Grid initial locations
    ind = sub2ind(size(tmask),ceil(init_x_all(S)),ceil(init_y_all(S)),floor(init_z_all(S)));
    age = nan(size(tmask));
    section = nan(size(tmask));
    volume = nan(size(tmask));
    age(ind) = final_age_all(S);
    section(ind) = final_section_all(S);
    volume(ind) = init_volume_all(S);
    % Chop them up and reduce precision
    age = single(age(region_limits(1,1):region_limits(1,2),region_limits(2,1):region_limits(2,2),:));
    section = single(section(region_limits(1,1):region_limits(1,2),region_limits(2,1):region_limits(2,2),:));
    volume = single(volume(region_limits(1,1):region_limits(1,2),region_limits(2,1):region_limits(2,2),:));
    
    % Tack them on
    age_all = cat(3,age_all,age);
    section_all = cat(3,section_all,section);
    volume_all = cat(3,volume_all,volume);
    
end
% Temporary arrays to remove young and unventilated particles
age_temp = age_all; age_temp(age_all<31536000 | section_all~=7)=NaN;
volume_temp = volume_all; volume_temp(age_all<31536000 | section_all~=7)=NaN;
% Volume-weighted mean
age_mean = nansum(volume_temp.*age_temp,3)./nansum(volume_temp,3);
% Find seeded locations
seeded_all = sum(isfinite(section_all),3);
% Set unventilated particles to -1
age_mean(sum(isfinite(age_temp),3)==0 & seeded_all~=0) = -1;
% Replace grids in which median age is less than 1 year with 0
age_mean(isnan(age_mean) & nanmedian(age_all,3)<31536000)=0; % All particles less than one year
clear age_temp volume_temp
% Pull them back together
age_meanT = nan(size(tmask(:,:,1)));
age_meanT(region_limits(1,1):region_limits(1,2),region_limits(2,1):region_limits(2,2))=age_mean;
age_mean = age_meanT;
clear age_all section_all volume_all

save('age_mean.mat','age_mean');


% CALCULATIONS
% Mean age in NADW
load('age_mean.mat');

% Mean latitude as a function of age
edges_age = 0:1:58;
latU = 70;

S = ariane_S({final_section_all,final_lat_all},{7,[-Inf latU]});
D = ariane_D(final_age_all/31536000,edges_age,S);
[Dwmean_lat,Dwstd_lat] = ariane_Dwmean(init_lat_all,init_volume_all,D);

% And split by density
S = ariane_S({final_section_all,final_lat_all,init_dens_all},{7,[-Inf latU],[27.7 27.8]});
D = ariane_D(final_age_all/31536000,edges_age,S);
[Dwmean_lat_dens78,Dwstd_lat_dens78] = ariane_Dwmean(init_lat_all,init_volume_all,D);

S = ariane_S({final_section_all,final_lat_all,init_dens_all},{7,[-Inf latU],[27.8 27.9]});
D = ariane_D(final_age_all/31536000,edges_age,S);
[Dwmean_lat_dens89,Dwstd_lat_dens89] = ariane_Dwmean(init_lat_all,init_volume_all,D);

S = ariane_S({final_section_all,final_lat_all,init_dens_all},{7,[-Inf latU],[27.9 28]});
D = ariane_D(final_age_all/31536000,edges_age,S);
[Dwmean_lat_dens90,Dwstd_lat_dens90] = ariane_Dwmean(init_lat_all,init_volume_all,D);

% FIGURES
close all
% Figure properties
gaph = 0.02;
lmap = 0.1; bmap = 0.5; wmap = 0.5; hmap = 0.4;
llat = lmap+wmap+gaph; blat = bmap; wlat = wmap/2; hlat = hmap;
posmap = [lmap bmap wmap hmap];
poslat = [llat blat wlat hlat];
clims = [-1 58]; nc = 58; cm = cmocean('matter',nc);
cm = [unventcolor; cm];
cmdens = cmocean('haline',3);
lw = 1.5;
 
figure('units','inches');
pos = get(gcf,'position');
set(gcf,'pos',[pos(1) pos(2) 8.3 11.7]);

axes('position',posmap);
worldmap(latlimB,lonlimB);
pcolorm(double(gphit),double(glamt),round(age_mean)/31536000); shading flat
framem('FlineWidth',2,'FEdgeColor','black');
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none');
colormap(cm); caxis(clims);
t = colorbar;
t.Limits = [0 58]; t.Label.String = 'Age (years)'; t.Label.Interpreter = 'latex';
t.Label.FontSize = axesfs;

title('(a) Age  of NADW (depth-mean)','fontsize',titlefs)

axes('position',poslat);
patch([centre(edges_age) fliplr(centre(edges_age))],...
   [Dwmean_lat+Dwstd_lat; flipud(Dwmean_lat-Dwstd_lat)],0.8*[1 1 1],'edgecolor','none');
hold on;
plot(centre(edges_age),Dwmean_lat,'w.');
h78 = plot(centre(edges_age),Dwmean_lat_dens78,'-','color',cmdens(1,:),'linewidth',lw);
h89 = plot(centre(edges_age),Dwmean_lat_dens89,'-','color',cmdens(2,:),'linewidth',lw);
h90 = plot(centre(edges_age),Dwmean_lat_dens90,'-','color',cmdens(3,:),'linewidth',lw);
set(gca,'yaxislocation','r','ylim',[-10 70],'xlim',[edges_age(1) edges_age(end)]);
box on; grid on;
xlabel('Age (years)'); ylabel('Latitude ($^\circ N$)');
l = legend([h78,h89,h90],{'27.7 to 27.8','27.8 to 27.9','27.9 to 28'},'location','southwest');
set(l,'box','off','interpreter','latex');
text(10,5,'$\gamma^N$:'); 
title('(b) Mean latitude as a function of age','fontsize',titlefs);

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_pathway'],'-png',saveres);
end

%% DIAPYCNAL MOTION
% CALCULATIONS
% Density difference as a function of age
edges_densdiff = -0.2:0.01:0.2; edges_age = 0:58;
Dsum_dage = nan(length(edges_densdiff)-1,length(edges_age)-1,length(times));
for t=1:length(times);
    disp(times(t));
    S = ariane_S({init_t_all,final_section_all, final_age_all},{times(t),7,[31536000 Inf]});
    D = ariane_D2(final_dens_all-init_dens_all,final_age_all/31536000,edges_densdiff,edges_age,S);
    [~,Dsum_age(:,:,t),~,~,~] = ariane_Dop(init_volume_all,D);
end
save('Dsum_dage.mat','Dsum_dage');

load('Dsum_dage.mat');
Dnum = 58:-1:1;
Dsummean_dage = nansum(Dsum_dage,3)./repmat(Dnum,[length(edges_densdiff)-1 1]);
Dnorm = Dsummean_dage./repmat(nansum(Dsummean_dage,1),[size(Dsummean_dage,1) 1]);

% Joint histogram
edges_dens = 27.7:0.01:28;
S = ariane_S({init_t,final_section, final_age},{4217,7,[31536000 Inf]});
D = ariane_D2(init_dens,final_dens,edges_dens,edges_dens,S);
[~,Dsum_d,~,~,~] = ariane_Dop(init_volume,D);

% Figure properties
close all;

clims_d = [0 1.5]; nc = 50; cm = cmocean('-gray',nc);
clims_dage = [0 0.15];
gapv = 0.1; gaph = 0.02;
lb = 0.2; bb = 0.1; wb = 0.6; hb = 0.4;
lcb = lb+wb+gaph; bcb = bb+0.1; wcb = 0.03; hcb = 0.2;
lt = lb; bt = bb+hb+gapv; wt = wb; ht = wt/2;
lct = lt+wt+gaph; bct = bt+0.05; wct = wcb; hct = hcb;
posb = [lb bb wb hb];
poscb = [lcb bcb wcb hcb];
post = [lt bt wt ht];
posct = [lct bct wct hct];

figure('units','inches');
pos = get(gcf,'position');
set(gcf,'pos',[pos(1) pos(2) 4.2 8.4]);

axes('position',post);
imagesc(centre(edges_dens),centre(edges_dens),Dsum_d'*1e-14); flipim;
hold on;
plot([edges_dens(1) edges_dens(end)],[edges_dens(1) edges_dens(end)],'k--');
axis square;
caxis(clims_d); colormap(cm);% colorbar;
xlabel('$\gamma^N$ in September 2015 ($\gamma^N_{2015}$)'); ylabel('$\gamma^N$ at subduction ($\gamma^N_{sub}$)');
text(27.975,27.725,'water becomes denser over time','horizontalalignment','right','fontsize',axesfs);
text(27.725,27.975,'water becomes lighter over time','horizontalalignment','left','fontsize',axesfs);
title('(a) Density changes in NADW','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(clims_d(1),clims_d(2),nc);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',posct)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cm);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('Volume ($10^{14}m^3$)','fontsize',axesfs)

axes('position',posb);
imagesc(centre(edges_densdiff),centre(edges_age),Dnorm');
flipim;
hold on
plot([0 0],[0 58],'k--');
text(0.05,55,'lighter','horizontalalignment','left','fontsize',axesfs);
text(-0.05,55,'denser','horizontalalignment','right','fontsize',axesfs)
xlabel('$\Delta \gamma^N = \gamma^N_{sub} - \gamma^N_{2015}$'); ylabel('Age (years)');
caxis(clims_dage); colormap(cm);
title('(b) Density changes as a function of age','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(clims_dage(1),clims_dage(2),nc);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscb)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cm);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('Fraction','fontsize',axesfs)

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_diapycnal'],'-png',saveres);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GIN SEAS RE-ENTRAINMENT
% WILL RE-RUN EXPERIMENT IN GLOBAL ORCA 5D, SEEDED IN SEPTEMBER

% Load data from GIN seas forwards-in-time experiment
model = '5d';
experiment = 'quant_for_seedGINn1_t1259_everywhere_MLrefz8delsig0.01';
fileout = [rootdir 'experiments/' model '/' experiment '/ariane_positions_quantitative.nc'];
filein = [rootdir 'experiments/' model '/' experiment '/ariane_initial.nc'];
gridf = 'ORCA025.L75-GJM189_mesh_hgr.nc';

final_x = ncread(fileout,'final_x');
final_y = ncread(fileout,'final_y');
init_x = ncread(fileout,'init_x');
init_y = ncread(fileout,'init_y');
init_volume = ncread(filein,'init_volume');
init_dens = ncread(filein,'init_dens');
final_dens = ncread(fileout,'final_dens');
final_section = ncread(fileout,'final_section');
final_lat = ncread(fileout,'final_lat');
final_age = ncread(fileout,'final_age');

glamt = ncread([rootdir 'grids_ncmod/' gridf ],'glamt');
gphit = ncread([rootdir 'grids_ncmod/' gridf ],'gphit');
[nx,ny]=size(glamt);

% Figure properties
gaph = 0.04; gapv = 0.08;
llat = 0.1; blat = 0.5; wlat = 0.2; hlat = 0.35;
lmap = llat+wlat+gaph; bmap = blat+gapv/2; wmap = wlat*2; hmap = hlat;
hrho = hmap/2; lrho = llat; brho = blat-gapv-hrho; wrho = wmap+gaph+wlat;

lcbarmap = lmap+gaph; bcbarmap = blat; wcbarmap = wmap-2*gaph; hcbarmap = gapv/4; 

posmap = [lmap bmap wmap hmap];
poslat = [llat blat wlat hlat];
posrho = [lrho brho wrho hrho];

poscbarmap = [lcbarmap bcbarmap wcbarmap hcbarmap];

cmrho = cbrewer('qual','Paired',12);

% Calculations
% Re-entrainment
Smap = ariane_S({final_section,final_age,init_x,init_y},{7,[31536000/4 Inf],[300 Inf],[440 Inf]});
edgesX = 0:nx;
edgesY = 0:ny;
D = ariane_D2(final_x,final_y,edgesX,edgesY,Smap);
[~,Dsum_map,~,~] = ariane_Dop(init_volume,D);
% Fraction of volume at each latitude
S = ariane_S({final_age},{[31536000/4 Inf]});
edges_lat = latlimB(1):2:latlimB(2);
D = ariane_D(final_lat,edges_lat,S);
[~,Dsum_lat,~,~] = ariane_Dop(init_volume,D);
Dsum_frac = Dsum_lat/sum(init_volume(S)); % Build fraction functionality into ariane_Dop
Dsum_frac(isnan(Dsum_frac))=0;
% Density shift
S = ariane_S({final_age,final_lat,final_section,init_x,init_y},{[31536000/4 Inf],[-Inf 67],7,[300 Inf],[440 Inf]});
edges_rho = 27.4:0.05:28.6;
Dinit = ariane_D(init_dens,edges_rho,S);
[~,Dinitsum,~,~] = ariane_Dop(init_volume, Dinit);
Dfinal = ariane_D(final_dens,edges_rho,S);
[~,Dfinalsum,~,~] = ariane_Dop(init_volume, Dfinal);

% Figure
close all
figure('units','inches');
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);

% Latitude
axes('position',poslat);
barh(centre(edges_lat),cumsum(Dsum_frac),1,'facecolor',cmrho(1,:));
ylim(latlimB);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
ylabel('Latitude ($^\circ N$)'); xlabel('Fraction');
grid on;
title('(a) Final latitude','fontsize',titlefs)

% Re-entrainment
axes('position',posmap);
worldmap(latlimN,lonlimN);
pcolorm(double(gphit),double(glamt),Dsum_map*1e-11);
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
l = legend('location','northwest','Initial density','Density at re-entrainment');
set(l,'interpreter','latex','box','off');

if strcmp(savefigs,'y')
    export_fig([savedir 'fig_GIN'],'-png',saveres);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MET OFFICE NADW SELECTION

% Load data from met office everywhere experiment
model = 'metoffice';
experiment = 'quant_back_seedNAn1_t4510_everywhere_MLrefz8delsig0.01';
fileout = [rootdir 'experiments/' model '/' experiment '/ariane_positions_quantitative.nc'];
filein = [rootdir 'experiments/' model '/' experiment '/ariane_initial.nc'];
gridf = 'anudlo_1m_mesh_hgr.nc';
load([rootdir 'time/time_' model '.mat'],'time');

init_dens = ncread(filein,'init_dens');
final_dens = ncread(fileout,'final_dens');
init_temp = ncread(fileout,'init_temp');
init_salt = ncread(fileout,'init_salt');
init_t = ncread(fileout,'init_t');
init_x = ncread(fileout,'init_x');
init_y = ncread(fileout,'init_y');
init_z = ncread(fileout,'init_z');
final_section = ncread(fileout,'final_section');
final_age = ncread(fileout,'final_age');
final_x = ncread(fileout,'final_x');
final_y = ncread(fileout,'final_y');
init_volume = ncread(filein,'init_volume');

% Grid details
tmask = ncread([rootdir 'grids_ncmod/' gridf ],'tmask');
glamt = ncread([rootdir 'grids_ncmod/' gridf ],'glamt');
gphit = ncread([rootdir 'grids_ncmod/' gridf ],'gphit');

% Density at final timestep
vosigntr = ncread([rootdir 'data_link/metoffice/anudlo_1m_4510_grid_R.nc'],'vosigntr');

% Determine YEAR of ventilation
[final_year,final_agereal] = calc_yr_ariane(final_age,init_t,model,time);
% Determine MONTH of ventilation
mtdec = (final_year-floor(final_year)); % Months as a decimal fraction of 1
        
% Put into discrete months
fm = 10; % fraction of month to divide into
final_month = nan(size(mtdec));
for mt = 1:12*fm;
    final_month(mtdec>=(mt-1)/(12*fm) & mtdec < mt/(12*fm))=mt/fm;
end;

% MONTH OF VENTILATION
close all
S = ariane_S({final_section,final_age, init_dens},{5,[31536000 Inf],[27.6 28]});

% Figure properties
lmap = 0.1; bmap = 0.5; wmap = 0.5; hmap = 0.4;
gapv = 0.03;
hhist = 0.1; lhist = lmap;  bhist = bmap-gapv-hhist; whist = wmap;
hcbar = 0.01; lcbar = lmap; bcbar = bhist-gapv/2-hcbar; wcbar = wmap;
posmap = [lmap bmap wmap hmap];
poshist = [lhist bhist whist hhist];
poscbar = [lcbar bcbar wcbar hcbar];
clims = [0 12]; nc = 12; cm = cbrewer('qual','Paired',nc);
cm = [unventcolor; cm];
 
figure('units','inches');
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

title('(a) Month of subduction in NADW (depth-median)','fontsize',titlefs)

% Histogram
edges = 0+1/(fm*2):1/fm:12+1/(fm*2);
D = ariane_D(final_month,edges,S);
[~,Dsum,~,~] = ariane_Dop(init_volume,D);
cmc = cbrewer('qual','Paired',nc);
axes('position',poshist);
for mt = 1:length(Dsum);
    fc = cmc(ceil(mt/fm),:);
    bar(edges(mt),Dsum(mt)*1e-15,1/fm,'facecolor',fc,'edgecolor','none');
    hold on
end
set(gca,'xticklabels',[]);
ylabel('Volume ($10^{15}\, m^3$)');
title('(b) Volume of NADW subducted in each month','fontsize',titlefs)

% Colour bar
yc = 0.5:12.5;
xc = 0:2;
[~,zc] = meshgrid(xc,yc);
ax = axes('position',poscbar);
contourf(yc,xc,zc',yc,'linestyle','none'); caxis(clims); colormap(ax,cbrewer('qual','Paired',12));
set(gca,'ytick',[],'xtick',1:12,'ticklength',[0 0],...
    'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
xlabel('Month of subduction')

% DENSITY DISTRIBUTION
edges_dens = 26:0.02:29;
S = ariane_S({final_section},{[-Inf Inf]}); % Capture all
Dall = ariane_D(init_dens,edges_dens,S);
[~,Dsumall,~,~] = ariane_Dop(init_volume,Dall);

S = ariane_S({final_section, final_age},{5, [31536000 Inf]});
D = ariane_D(init_dens,edges_dens,S);
[~,Dsum,~,~] = ariane_Dop(init_volume,D);

close all;
figure;
bar(centre(edges_dens),Dsumall,1,'facecolor',0.4*[1 1 1],'edgecolor','none');
hold on;
bar(centre(edges_dens),Dsum,1,'facecolor',0.8*[1 1 1],'edgecolor','none');

% SUBDUCTION LOCATION
S = ariane_S({final_section, final_age, init_dens},{5, [31104000 50*31104000],[27 27.5]});

minX = floor(min(final_x(S))); maxX = ceil(max(final_x(S)));
minY = floor(min(final_y(S))); maxY = ceil(max(final_y(S)));
edgesX = minX:1:maxX;
edgesY = minY:1:maxY;
D = ariane_D2(final_x,final_y,edgesX,edgesY,S);
[Dind,Dsum,~,~] = ariane_Dop(init_volume,D);
Dsum(Dsum==0)=NaN;

TrangeX = minX+1:maxX; TrangeY = minY+1:maxY;

close all;
figure
%%% Location of subduction
worldmap(latlimB,lonlimB);
pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),Dsum*1e-13); shading flat
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmsub); caxis(climssub)

% Grid initial locations
ind = sub2ind(size(tmask),ceil(init_x),ceil(init_y),floor(init_z));
volume = nan(size(tmask));
volume(ind(S)) = init_volume(S);
volume_sum = nansum(volume,3);

figure
%%% Ventilated location
worldmap(latlimB,lonlimB);
pcolorm(double(gphit),double(glamt),volume_sum); shading flat
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmsub);% caxis(climssub)

%%% DENSITY CHANGES
t1 = 1; t2 = Inf;
edges_dens = 27:0.02:28;
S = ariane_S({final_section, final_age},{5,[t1*31104000 t2*31104000]});
D = ariane_D2(init_dens,final_dens,edges_dens,edges_dens,S);
[Dind,Dsum,~,~] = ariane_Dop(init_volume,D);

%close all;
figure;
pcolor(centre(edges_dens),centre(edges_dens),Dsum./sum(init_volume(S)));
shading flat;
axis square
title([num2str(t1) ' to ' num2str(t2)]);
caxis([0 1]*0.01)
hold on
plot([edges_dens(1) edges_dens(end)], [edges_dens(1) edges_dens(end)], 'w--');

% Density difference as a function of age
close all 
edges_densdiff = -0.5:0.02:0.5; edges_age = 0:1/12:1;%370; 
S = ariane_S({final_section, final_dens},{5,[27 27.7]});
D = ariane_D2(final_dens-init_dens,final_age/31104000,edges_densdiff,edges_age,S);
[~,Dsum,~,~] = ariane_Dop(init_volume,D);

imagesc(centre(edges_densdiff),centre(edges_age),Dsum'); shading flat;

% Density stratification
close all;
x = 1170;
figure;
contourf(squeeze(vosigntr(x,300:900,:))',22:0.1:28.5,'edgecolor','none');
hold on;
contour(squeeze(vosigntr(x,300:900,:))',27.6:0.1:28,'k-');
set(gca,'ydir','r')

close all
universal
savefigs='y';
if exist('overrulesave')==1;
	savefigs=overrulesave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL NADW SELECTION

if donenadwcalcs=='n';
	% Load results of everywhere ariane experiment
	init_t_e = ncread(fileout_e,'init_t');
	final_section_e = ncread(fileout_e,'final_section');
	final_lat_e = ncread(fileout_e,'final_lat');
	init_dens_e = ncread(fileout_e,'init_dens');
	init_volume_e = ncread(filein_e,'init_volume');
	
	rr = 27:0.02:max(init_dens_e);
	rrC = 0.5*(rr(1:end-1)+rr(2:end));
	S = ariane_S({init_t_e},{4217});
	Sv = ariane_S({init_t_e,final_section_e},{4217,7});	
	
	D = ariane_D(init_dens_e,rr,S);
	Dv = ariane_D(init_dens_e,rr,Sv);
	
	[~,Dsum,~,~] = ariane_Dop(init_volume_e,D);
	[~,Dvsum,~,~] = ariane_Dop(init_volume_e,Dv);

	% Also load density
	vosigntr = calc_sigmantr(ncread(fileden,'votemper'),ncread(fileden,'vosaline'));;
	vosigntr(tmask==0)=NaN;
	vosigntr_xmean = squeeze(nanmean(vosigntr,1));

	%[nx,ny,nz]=size(vosigntr);
	% Calculate the AMOC
	disp('calculating AMOC')
	AMOC = nan(ny,nz);
	AMOCdens = nan(ny,length(rr)-1);
        for j = 1:ny;
            vomecrty = squeeze(ncread(filevel,'vomecrty',[1 j 1 1],[Inf 1 Inf Inf]));
            e1v_j = squeeze(e1v(:,j,:)); e3v_j = squeeze(e3v(:,j,:));
	    tmask_j = squeeze(tmask(:,j,:));
	    V = vomecrty.*e1v_j.*e3v_j.*tmask_j;
	    
	    vosigntr_j = squeeze(vosigntr(:,j,:));
	    Vdens = nan(length(rr)-1,1);
	    for i=1:length(rr)-1;
		    Vtemp = V(vosigntr_j>rr(i) & vosigntr_j<=rr(i+1));
	 	    Vdens(i) = nansum(Vtemp(:));
	    end

            AMOC(j,:)=cumsum(nansum(V,1),'reverse');
            AMOCdens(j,:)=cumsum(Vdens,'reverse');
	end

	resetcalcs
	donenadwcalcs='y';	
end

% Figure properties
gaph=0.05; gapv=0.05;
lc=0.1; bc=0.1; wc=0.5; hc=0.15;
lb=lc; bb=bc+hc+gapv; wb=wc; hb=hc;
la=lb; ba=bb+hb+gapv; wa=wb; ha=hb*2;
wcbar = 0.02;

posssd=[la,ba,wa,ha];
posven=[lb,bb,wb,hb];
posmoc=[lc,bc,wc,hc];

poscbarssd = [la+wa-gaph,ba+ha/2,wcbar,ha/3];
poscbarmoc = [lc+wc+gaph/2,bc,wcbar,hc];

% Overturning circulation in density space
f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 8.3]);

axmoc = axes('position',posmoc);
title('Overturning circulation in density space (late winter mean)','fontsize',titlefs)
x  = gphit(1000,:);
x(end)=x(end-1)+abs(diff([x(end),x(end-1)]));
contourf(x,rrC,-AMOCdens'*1e-6,ncmoc,'LineColor','none');
hold on;
plot(x,27.7*ones(length(x),1),'w--')
plot(x,28*ones(length(x),1),'w--')
set(gca,'ydir','r')
xlim([20 70]);
ylim([27 28.1])
xlabel('Latitude ($^\circ N$)')
ylabel('Density ($kgm^{-3}$)')
caxis(climsmoc);
colormap(gca,cmmoc);

% Colorbar
xc = 0:2;
yc = linspace(climsmoc(1),climsmoc(2),ncmoc);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarmoc)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cmmoc);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('$Sv$','fontsize',axesfs)

if strcmp(savefigs,'y')
    	disp('SAVING FIGURE')
	%export_fig([savedir 'fig_nadw'],'-png',saveres);
	%saveas(f,[savedir 'fig_nadw'],'png');
	print(f, [savedir 'fig_amoc_dens'], '-dpng',saveres);
end
close all

f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);

% Surface density
axes('position',posssd)
worldmap(latlimN,lonlimN);
pcolorm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),...
	vosigntr_wmean(TrangeX,TrangeY,8)); shading flat
contourm(double(gphit(TrangeX,TrangeY)),double(glamt(TrangeX,TrangeY)),...
	vosigntr_wmean(TrangeX,TrangeY,8),[27.7, 27.8, 27.9, 28],'color','k')
geoshow('landareas.shp','facecolor',landcolor,'edgecolor','none')
framem('FlineWidth',2,'FEdgeColor','black');
colormap(gca,cmssd); caxis(climsssd)

title('(a) Neutral density (late-winter mean) at 10m','fontsize',titlefs)

% Colorbar
xc = 0:2;
yc = linspace(climsssd(1),climsssd(2),ncssd);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarssd)
contourf(xc,yc,zc,yc,'LineColor','none');
hold on
for i=27.7:0.1:28;
	plot(xc,[i,i,i],'k-')
end
colormap(gca,cmssd);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('$kg\,m^{-3}$','fontsize',axesfs)

% Ventilated distribution

axes('position',posven)
bar(rrC,Dsum/1e15,1,'facecolor',0.8*[1 1 1],'edgecolor','none');
hold on;
bar(rrC,Dvsum/1e15,1,'facecolor',0.6*[1 1 1],'edgecolor','none');
ylim([0 16]);
xlim([27 28.4]);
xlabel('$\gamma^n$ ($kgm^{-3}$)')
ylabel('Volume ($10^{15}\,m^3$)')
l=legend('Total','Ventilated','location','northwest');
set(l,'interpreter','latex','box','off')

title('(b) Density distribution of Atlantic basin and ventilated fraction','fontsize',titlefs)

% Overturning circulation
axmoc = axes('position',posmoc)
title('(c) Overturning circulation (late winter mean)','fontsize',titlefs)
x  = gphit(1000,:);
x(end)=x(end-1)+abs(diff([x(end),x(end-1)]));
contourf(x,gdept,-AMOC'*1e-6,ncmoc,'LineColor','none');
hold on;
contour(x,gdept,vosigntr_xmean',[27.7 28])
set(gca,'ydir','r')
xlim([20 70]);
xlabel('Latitude ($^\circ N$)')
ylabel('Depth ($m$)')
caxis(climsmoc);
colormap(gca,cmmoc);

axes('position',posmoc)
b = vosigntr_xmean;
b(isnan(b))=-1;
b(b>0)=NaN;
pcolor(x,gdept,b'); shading flat;
set(gca,'ydir','r');
xlim([20 70]);
caxis([-2 1]);
colormap(gca,copper(4))
set(gca,'visible','off')

title(axmoc,'(c) Overturning circulation (late winter mean)','fontsize',titlefs)
% Colorbar
xc = 0:2;
yc = linspace(climsmoc(1),climsmoc(2),ncmoc);
[~,ycg]=meshgrid(xc,yc);
zc=ycg;

axes('position',poscbarmoc)
contourf(xc,yc,zc,yc,'LineColor','none');
colormap(gca,cmmoc);
set(gca,'xtick',[],'yaxislocation','right');
ylabel('$Sv$','fontsize',axesfs)

if strcmp(savefigs,'y')
    	disp('SAVING FIGURE')
	%export_fig([savedir 'fig_nadw'],'-png',saveres);
	%saveas(f,[savedir 'fig_nadw'],'png');
	print(f, [savedir 'fig_nadw'], '-dpng',saveres);
end


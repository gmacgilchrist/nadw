close all
universal
savefigs='y';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MLD ALONG PSI, SUBDUCTION VOLUME ALONG PSI
if donemldonpsicalcs=='n'
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
	save('matfiles/D2_fXfY_sum_volume.mat','D','Dsum');
	else
	load('matfiles/D2_fXfY_sum_volume.mat');
	end
	
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

	resetcalcs
	donemldonpsicalcs='y'
end
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

f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 8.3 11.7]);
axes('position',posmap);
axmap = worldmap(latlimLS,lonlimLS);
axmld = axes('position',posmld, 'ydir','r'); hold(axmld,'on')
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
    if c == 53;
    	is = [50,200];
	for d = is;
		axes(axmap)
		textm(double(clat{c}(ind(d))),double(clon{c}(ind(d))),num2str(round(dist(d))));
	end
    end
    plot(axmld,dist,cmld{c}(ind),'color',cmpsi(40+clev_subp{c},:));
    plot(axsub,dist,cDsum{c}(ind)*1e-13,'color',cmpsi(40+clev_subp{c},:));
end
set(axmld,'ydir','r');
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
xlabel(axsub,'Distance along streamline ($km$)'); ylabel(axsub,'Volume ($10^{13}m^3$)');
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
    %export_fig([savedir 'fig_mld-on-psi'],'-png',saveres);
%	saveas(f,[savedir 'fig_mld-on-psi'],'png');
	print(f, [savedir 'fig_mld-on-psi'], '-dpng',saveres);
end


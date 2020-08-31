close all
universal
savefigs='y';
if exist('overrulesave')==1;
        savefigs=overrulesave;
end

%% DIAPYCNAL MOTION
% CALCULATIONS
if donediacalcs=='n';
	edges_densdiff = -0.2:0.01:0.2; edges_age = 0:58;
	if ~exist('matfiles/Dsum_dage.mat','file')
		% Density difference as a function of age
		Dsum_dage = nan(length(edges_densdiff)-1,length(edges_age)-1,length(times));
		for t=1:length(times);
	    		disp(times(t));
	    		S = ariane_S({init_t_all,final_section_all, final_age_all},{times(t),7,[31536000 Inf]});
	    		D = ariane_D2(final_dens_all-init_dens_all,final_age_all/31536000,edges_densdiff,edges_age,S);
	    		[~,Dsum_age(:,:,t),~,~,~] = ariane_Dop(init_volume_all,D);
		end
		save('Dsum_dage.mat','Dsum_dage');
	else
		load('matfiles/Dsum_dage.mat');
	end
	Dnum = 58:-1:1;
	Dsummean_dage = nansum(Dsum_dage,3)./repmat(Dnum,[length(edges_densdiff)-1 1]);
	Dnorm = Dsummean_dage./repmat(nansum(Dsummean_dage,1),[size(Dsummean_dage,1) 1]);
	
	% Joint histogram
	edges_dens = 27.7:0.01:28;
	S = ariane_S({init_t,final_section, final_age},{4217,7,[31536000 Inf]});
	D = ariane_D2(init_dens,final_dens,edges_dens,edges_dens,S);
	[~,Dsum_d,~,~,~] = ariane_Dop(init_volume,D);
	
	resetcalcs
	donediacalcs='y'
end

% FIGURE
% Figure properties

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

f = figure('units','inches');
if strcmp(savefigs,'y');
	set(f,'visible','off');
end
pos = get(gcf,'position');
set(gcf,'pos',[pos(1) pos(2) 4.2 8.4]);

axes('position',post);
imagesc(centre(edges_dens),centre(edges_dens),Dsum_d'*1e-14); flipim;
hold on;
plot([edges_dens(1) edges_dens(end)],[edges_dens(1) edges_dens(end)],'k--');
axis square;
caxis(clims_d); colormap(cm);% colorbar;
xlabel('$\gamma^N$ in September 2015 ($\gamma^N_{2015}$)'); ylabel('$\gamma^N$ at subduction ($\gamma^N_{sub}$)');
text(27.975,27.725,'water gets denser over time','horizontalalignment','right','fontsize',axesfs);
text(27.725,27.975,'water gets lighter over time','horizontalalignment','left','fontsize',axesfs);
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
%flipim;
hold on
plot([0 0],[0 58],'k--');
text(0.05,5,'gets lighter','horizontalalignment','left','fontsize',axesfs);
text(-0.05,5,'gets denser','horizontalalignment','right','fontsize',axesfs)
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
    %export_fig([savedir 'fig_diapycnal'],'-png',saveres);
%	saveas(f,[savedir 'fig_diapycnal'],'png');
	print(f, [savedir 'fig_diapycnal'], '-dpng',saveres);
end

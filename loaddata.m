% Load all of the data
% First check if the variables exist already

if doneloaddata=='n'
	disp('LOADING DATA')
	% model data
	load([rootdir 'time/time_' model '.mat'],'time');
	tmask = single(ncread([rootdir 'grids_ncmod/' gridf ],'tmask'));
	gphit = ncread([rootdir 'grids_ncmod/' gridf ],'gphit');
	glamt = ncread([rootdir 'grids_ncmod/' gridf ],'glamt');
	e1t = ncread([rootdir 'grids_ncmod/' gridf ],'e1t');
	e2t = ncread([rootdir 'grids_ncmod/' gridf ],'e2t');
	e3t = ncread([rootdir 'grids_ncmod/' gridf ],'e3t');
	e1v = repmat(ncread([rootdir 'grids_ncmod/' gridf ],'e1v'),[1 1 75]);
	e3v = ncread([rootdir 'grids_ncmod/' gridf ],'e3t');
	V = repmat(e1t,[1 1 75]).*repmat(e2t,[1 1 75]).*e3t;
	% model data direct from model
	modeldir = '/home/ocean_shared_data1/DRAKKAR/ORCA025.L75-GJM189-S/';
	meanfile = 'MarchApril_means_1975-2015/ORCA025.L75-GJM189_Mar-Apr-1975-2015-mean_';
	somxl010_wmean = ncread([modeldir meanfile 'gridT.nc'],'somxl010');
	vovecrtz_wmean = ncread([modeldir meanfile 'gridW.nc'],'vovecrtz');
        vosaline_wmean = ncread([modeldir meanfile 'gridT.nc'],'vosaline');
        votemper_wmean = ncread([modeldir meanfile 'gridT.nc'],'votemper');
	sobarstf_wmean = ncread([modeldir meanfile 'psi.nc'],'sobarstf');
	gridfile = 'GRID/ORCA025.L75-GJM189_mesh_zgr.nc';
	mbathy = ncread([modeldir gridfile],'mbathy');
	gdept = ncread([modeldir gridfile],'gdept_0');
	% neutral density of gridded data
	vosigntr_wmean = calc_sigmantr(votemper_wmean,vosaline_wmean);

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
	
	doneloaddata='y'
end

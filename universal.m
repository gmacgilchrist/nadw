% Anything that is common to all figures
% e.g. Figure specs and file paths

savedir = '/home/ocean_personal_data/graemem/frontend/paperfigures/2018_nadw/';
saveres='-r450';

%%
% file locations
rootdir = '/home/ocean2/graemem/ariane/';
model = 'orca025_global_5d';
experiment = 'quant_back_seedNAn1_t3560-sep-4217_sign27.7-28_MLrefz8delsig0.01';
fileout = [rootdir 'experiments/' model '/' experiment '/ariane_positions_quantitative.nc'];
filein = [rootdir 'experiments/' model '/' experiment '/ariane_initial.nc'];
gridf = 'ORCA025.L75-GRD88_mesh_hgr.nc';
gridzf = 'ORCA025.L75-GRD88_mesh_zgr.nc';
%'everywhere' experiment
experiment_e = 'quant_back_seedNAn1_t4217_everywhere_MLrefz8delsig0.01';
fileout_e = [rootdir 'experiments/' model '/' experiment_e '/ariane_positions_quantitative.nc'];
filein_e = [rootdir 'experiments/' model '/' experiment_e '/ariane_initial.nc'];
% GIN experiment
model_5d = '5d';
experiment_g = 'quant_for_seedGINn1_t1259_everywhere_MLrefz8delsig0.01';
fileout_g = [rootdir 'experiments/' model_5d '/' experiment_g '/ariane_positions_quantitative.nc'];
filein_g = [rootdir 'experiments/' model_5d '/' experiment_g '/ariane_initial.nc'];
gridf_g = 'ORCA025.L75-GJM189_mesh_hgr.nc';
% Velocity for AMOC calculation
model = 'orca025_global_5d';
config = 'ORCA025.L75-GJM189';
filevel = [rootdir 'data_link/' model '/' config '_Mar-Apr-1975-2015-mean_gridV.nc'];
fileden = [rootdir 'data_link/' model '/' config '_Mar-Apr-1975-2015-mean_gridT.nc'];
filew = [rootdir 'data_link/' model '/' config '_Mar-Apr-1975-2015-mean_gridW.nc'];
% qualitiative experiment - diapycnal
model = 'orca025_global_5d';
experiment_d78 = 'qual_back_seedNAn1_t3560-sep-4217_sign27.7-28_MLrefz8delsig0.01_subbin-t4217-sec7-sig27.7-27.8';
experiment_d89 = 'qual_back_seedNAn1_t3560-sep-4217_sign27.7-28_MLrefz8delsig0.01_subbin-t4217-sec7-sig27.8-27.9';
experiment_d90 = 'qual_back_seedNAn1_t3560-sep-4217_sign27.7-28_MLrefz8delsig0.01_subbin-t4217-sec7-sig27.9-28';
filemask_d78 = [rootdir 'experiments/' model '/' experiment_d78 '/matfiles/traj_truncated_quant_MLrefz8delsig0.01.mat'];
filemask_d89 = [rootdir 'experiments/' model '/' experiment_d89 '/matfiles/traj_truncated_quant_MLrefz8delsig0.01.mat'];
filemask_d90 = [rootdir 'experiments/' model '/' experiment_d90 '/matfiles/traj_truncated_quant_MLrefz8delsig0.01.mat'];
filemap_d78 = [rootdir 'experiments/' model '/' experiment_d78 '/matfiles/criter1_MLrefz8delsig0.01/D2_iUjV_Dmean_drdt.mat'];
filemap_d89 = [rootdir 'experiments/' model '/' experiment_d89 '/matfiles/criter1_MLrefz8delsig0.01/D2_iUjV_Dmean_drdt.mat'];
filemap_d90 = [rootdir 'experiments/' model '/' experiment_d90 '/matfiles/criter1_MLrefz8delsig0.01/D2_iUjV_Dmean_drdt.mat'];

filedia = [rootdir 'experiments/' model '/qual_exps_combinedcalcs/D2_iUjV_Dmean_drdt.mat'];
filedia_onlyD = [rootdir 'experiments/' model '/qual_exps_combinedcalcs/D2_iUjV_Dmean_drdt_onlyD.mat'];
filekz = [rootdir 'experiments/' model '/qual_exps_combinedcalcs/D2_kzdens_Dsum_volume.mat'];
%%
% Figure specs
% Mask colors
landcolor = 0.6*[1 1 1];
unventcolor = 0.4*[1 1 1];

% Basin-wide map limits
latlimB = [-30 80]; lonlimB = [-100 20];
% Northern region map limits
latlimN = [45 75]; lonlimN = [-65 10];
% Northern region + GIN  map limits
latlimNG = [45 80]; lonlimNG = [-65 20];
% Labrador Sea map limits
latlimLS = [48 65]; lonlimLS = [-62 -42];
% Northern hemishpere basin
latlimNH = [-30 70]; lonlimNH = [-100 0];
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
% Mixed layer density
climsssd = [25,28.5]; ncssd = 50; cmssd = cmocean('haline',ncssd);
% AMOC
climsmoc = [-20 20]; ncmoc = 20; cmmoc = cmocean('balance',ncmoc);
% kz
climskz = [-6 1]; nckz = 8; cmkz = cbrewer('seq','Blues',nckz);

loaddata
globalcalculations


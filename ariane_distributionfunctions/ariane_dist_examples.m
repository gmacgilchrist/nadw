% Example codes for evaluating distributions

% ORDER OF OPERATION
% 1. Load ariane output
% 2. Define broad subsetting categories, e.g. only ventilated particles
% 3. Split this subset into a distribution based on certain
% characteristics, e.g. particle age
% 4. Perform operation on this distribution, e.g. sum volumes in each bin
%

clear; close all;
rootdir = '/home/ocean2/graemem/ariane/experiments/';
model = 'orca025_global_5d';
experiment = 'quant_back_seedNAn1_t3560-sep-4217_sign27.7-28_MLrefz8delsig0.01';
fileout = [rootdir model '/' experiment '/ariane_positions_quantitative.nc'];
fileinit = [rootdir model '/' experiment '/ariane_initial.nc'];
final_section = ncread(fileout,'final_section');
init_t = ncread(fileout,'init_t');
final_age = ncread(fileout,'final_age');
init_temp = ncread(fileout,'init_temp');
init_volume = ncread(fileinit,'init_volume');
final_x = ncread(fileout,'final_x');
final_y = ncread(fileout,'final_y');

%% VOLUMETRIC AGE DISTRIBUTION OF VENTILATED PARTICLES AT t=3980
% D_a(V(a>1yr,t=3980));

subchars={final_section,final_age,init_t};
values={7,[31536000/2 Inf],4217};
S=ariane_S(subchars,values);

distchar = final_age/31536000;
edges = 0:1:55;

D = ariane_D(distchar, edges, S);

% Age distribution summed by volume
sumchar = init_volume;
[~,Dsum,~,~] = ariane_Dop(sumchar,D);
figure; bar(centre(D.edges),Dsum)

%% AVERAGE (OVER INITIAL TIME) VOLUMETRIC AGE DISTRIBUTION OF VENTILATED PARTICLES
times = unique(init_t);

subchars={final_section,final_age};
values={7,[0 Inf]};
S=ariane_S(subchars,values);
D = ariane_D2(final_age/31536000,init_t,0:55,times-0.5,S);
Dsum = ariane_Dop(init_volume,D,'sum');

meanDsum = sum(Dsum,2)./sum(isfinite(Dsum),2);
figure; bar(meanDsum);

%% VOLUMETRIC DISTRIBUTION OF VENTILATION LOCATION
close all

subchars={final_section,final_age,init_t};
values={7,[31536000/2 Inf],4217};
S=ariane_S(subchars,values);

D = ariane_D2(final_x,final_y,floor(min(final_x)):ceil(max(final_x)),floor(min(final_y)):ceil(max(final_y)),S);
[Dind,Dsum,Dmean,Dmedian] = ariane_Dop(init_volume,D);
figure; imagesc(Dsum'); flipim;

% Do it through indexing
v = init_volume(S);
[nx,ny]=size(Dind);
for i=1:nx
    for j=1:ny
        vsum(i,j)=sum(v(Dind{i,j}(:)));
    end
end
figure; imagesc(vsum'); flipim;
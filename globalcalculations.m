% Do some preliminary global calculations
% Don't redo if they've already been done

if doneglobalcalcs=='n'
	disp('PERFORMING GLOBAL CALCULATIONS')
	
	[nx,ny,nz]=size(e3t);
	% Mask
	tmask(:,944:end,:)=0; % Very northern edge of domain
	tmask(1234:1244,924:926,:)=0; % Fiddly bit in the NorthEast
	tmask(1:893,839:end,:)=0; % CAA
	tmask(893:908,866:943,:)=0; % CAA
	tmask(1:862,754:838,:)=0; % Rest of Hudson Bay
	tmask(971:975,942:943,:)=0; % Fiddly bit in Baffin Bay
	tmask(1228:end,1:627,:)=0; % Indian Ocean and Red Sea
	tmask(1:758,:,:)=0; % North Pacific
	tmask(1:919,1:525,:)=0; % South Pacific
	% Pacific side of Central America
	for i=759:839;
	   ind = find(tmask(i,526:end,:)<1,1,'first')+525;
	   tmask(i,1:ind,:)=0;
	end
	% Fiddly bits in Central America
	tmask(829,531:532,:)=0;
	tmask(823,530:531,:)=0;
	tmask(807,542:543,:)=0;
	tmask(1147:end,627:731,:)=0; % Western Med
	tmask(1126:end,650:662,:)=0; % Gibraltar
	tmask(:,1:351,:)=0; % South Atlantic
	% Great Lakes
	tmask(784:846,687:739,:)=0;

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
	%27.7-27.8
	ind = find(init_t == 4217 & final_section==7 & init_dens>27.7 & init_dens<27.8);
	fileid = fopen([rootdir 'experiments/' model '/qual' experiment(6:end) '_subbin-t4217-sec7-sig27.7-27.8/subset.txt'],'w');
	fprintf(fileid,'%i \n',ind);
	fclose(fileid);
	%27.8-27.9
	ind = find(init_t == 4217 & final_section==7 & init_dens>=27.8 & init_dens<27.9);
        fileid = fopen([rootdir 'experiments/' model '/qual' experiment(6:end) '_subbin-t4217-sec7-sig27.8-27.9/subset.txt'],'w');
        fprintf(fileid,'%i \n',ind);
        fclose(fileid);
	%27.9-28
	ind = find(init_t == 4217 & final_section==7 & init_dens>=27.9 & init_dens<=28);
        fileid = fopen([rootdir 'experiments/' model '/qual' experiment(6:end) '_subbin-t4217-sec7-sig27.9-28/subset.txt'],'w');
        fprintf(fileid,'%i \n',ind);
        fclose(fileid);

	% Distribution at subduction location
        S = ariane_S({init_t,final_section,final_age},{4217,7,[31536000 Inf]});
	minX = floor(min(final_x(S))); maxX = ceil(max(final_x(S)));
        minY = floor(min(final_y(S))); maxY = ceil(max(final_y(S)));
        % Range of tracer points on model grid (shifted right and up for
        % correspondence with model indexing)
        TrangeX = minX+1:maxX; TrangeY = minY+1:maxY;
	
	doneglobalcalcs='y'
end

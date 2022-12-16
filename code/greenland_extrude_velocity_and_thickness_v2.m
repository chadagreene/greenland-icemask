% This script takes extrapolates ITS_LIVE v2 velocity and BedMachine v5
% thickness beyond the current extents of Greenland Ice Sheet. 
% 
% Written by Chad A. Greene, July 2022, updated November 2022. 
% 
% This script takes about an hour to run on my laptop from 2019. 

%% Load data 

% Load ITS_LIVE velocity grid: 
fn = 'ITS_LIVE_velocity_120m_GRE_0000_v02.nc';
x = ncread(fn,'x'); 
y = ncread(fn,'y'); 
vx = permute(ncread(fn,'vx0'),[2 1]);
vy = permute(ncread(fn,'vy0'),[2 1]);
v_error = permute(ncread(fn,'v0_error'),[2 1]);
outlier_frac = permute(ncread(fn,'outlier_frac'),[2 1]);
v_amp_error = permute(ncread(fn,'v_amp_error'),[2 1]); % an extra metric of uncertainty 
[X,Y] = meshgrid(x,y); 

% Get rid of bad velocity data: 
dst = bwdist(isnan(vx));
bad = v_error>10 | outlier_frac>0.05 | hypot(vx,vy)>25e3 | v_amp_error>10 | isnan(v_error) | bwdist(isnan(vx))<5; 
bad(v_error>3 & dst<10) = true; 

vx(bad) = nan; 
vy(bad) = nan; 
v_error(bad) = nan; 

fn = '/Users/cgreene/Documents/data/velocity/its_live/GRE_G0120_0000.nc';
x_its = ncread(fn,'x');
y_its = ncread(fn,'y');
v_error_its = permute(ncread(fn,'v_err'),[2 1]);
chip_its = permute(ncread(fn,'chip_size_max'),[2 1]);
count_its = permute(ncread(fn,'count'),[2 1]);
vx_its = permute(ncread(fn,'vx'),[2 1]);
vy_its = permute(ncread(fn,'vy'),[2 1]);
dst = bwdist(vx_its==-32767);

bad = v_error_its==0 | v_error_its>20 | hypot(vx_its,vy_its)>25e3; 
bad(chip_its<1920 & v_error_its>3 & dst<10) = true; % slightly high error places near the edge of the data is mostly in fjords and mostly bad. 
bad(count_its<100) = true; 

v_error_its(bad) = nan; 
vx_its(bad) = nan; 
vy_its(bad) = nan; 

v_error_its_interp = interp2(x_its,y_its,v_error_its,X,Y); 
isn = isnan(v_error_its_interp); 
v_error_its_interp(isn) = interp2(x_its,y_its,v_error_its,X(isn),Y(isn),'nearest'); 

vy_its_interp = interp2(x_its,y_its,vy_its,X,Y); 
isn = isnan(vy_its_interp); 
vy_its_interp(isn) = interp2(x_its,y_its,vy_its,X(isn),Y(isn),'nearest'); 

vx_its_interp = interp2(x_its,y_its,vx_its,X,Y); 
isn = isnan(vx_its_interp); 
vx_its_interp(isn) = interp2(x_its,y_its,vx_its,X(isn),Y(isn),'nearest'); 


use_v1 = isnan(vx) | isnan(v_error) | v_error_its_interp<v_error; 
vx(use_v1) = vx_its_interp(use_v1); 
vy(use_v1) = vy_its_interp(use_v1); 
v_error(use_v1) = v_error_its_interp(use_v1); 

clear bad chip_its count_its dst fn isn outlier_frac use_v1 v v_error_its* vx_its* vy_its* x_its y_its v_amp_error

% v_source = 5*ones(size(vx),'uint8'); % 5 = extrapolated
% v_source(isfinite(vx)) = 1;          % 1 = observed

% Joughin's data: 
[tmp,xj,yj] = geoimread('greenland_vel_mosaic250_vx_v1.tif'); 
tmp(abs(tmp)>25e3) = nan; 
vxj = interp2(xj,yj,tmp,X,Y); 

tmp = imread('greenland_vel_mosaic250_vy_v1.tif');
tmp(abs(tmp)>25e3) = nan; 
vyj = interp2(xj,yj,tmp,X,Y); 

tmp = imread('greenland_vel_mosaic250_ex_v1.tif');
tmp(abs(tmp)>25e3) = nan; 
vx_errorj = interp2(xj,yj,tmp,X,Y); 

tmp = imread('greenland_vel_mosaic250_ey_v1.tif');
tmp(abs(tmp)>25e3) = nan; 
vy_errorj = interp2(xj,yj,tmp,X,Y); 

v_errorj = hypot(vx_errorj.*vxj./hypot(vxj,vyj),vy_errorj.*vyj./hypot(vxj,vyj));

weighti = 1./v_error.^2; 
weightj = 1./v_errorj.^2; 

vx_tmp = vx; 
vy_tmp = vy; 

vx = nan(size(vx)); 
vx(isfinite(vx_tmp) & isnan(vxj)) = vx_tmp(isfinite(vx_tmp) & isnan(vxj)); 
vy = nan(size(vx)); 
vy(isfinite(vy_tmp) & isnan(vxj)) = vy_tmp(isfinite(vy_tmp) & isnan(vxj)); 

vx(isnan(vx_tmp) & isfinite(vxj)) = vxj(isnan(vx_tmp) & isfinite(vxj)); 
vy(isnan(vy_tmp) & isfinite(vxj)) = vyj(isnan(vy_tmp) & isfinite(vxj)); 

both = isfinite(weighti) & isfinite(weightj) & weighti>0 & weightj>0;  
vx(both) = (vx_tmp(both).*weighti(both) + vxj(both).*weightj(both))./(weighti(both) + weightj(both)); 
vy(both) = (vy_tmp(both).*weighti(both) + vyj(both).*weightj(both))./(weighti(both) + weightj(both)); 

% Islands: 
bad = inpolygon_map(double(x),double(y),[-426191     -246815     -272921     -443034],[-964839     -829254     -752620     -860414]);
vx(bad) = nan; 
vy(bad) = nan; 

clear vx_tmp vy_tmp both weighti weightj bad xj yj vxj vyj vx_errorj vy_errorj tmp v_errorj v_error

% Bedmachine thickness: 
[th_bm,x_bm,y_bm] = bedmachine_data('thickness','greenland'); 
bed_bm = bedmachine_data('bed','greenland'); 
mask_bm = bedmachine_data('mask','greenland'); 
errbed_bm = bedmachine_data('errbed','greenland'); 
errbed_bm(mask_bm==0) = hypot(1.12*1.047*errbed_bm(mask_bm==0),32.7); % 32.7 m is the inversion_error value from thickness_inversion_uncertainty. 
th_bm_inv = -1.12*1.047*bed_bm; % 1.12 converts bed to thickness at hydrostatic equilibrium and 1.047 is bc the median bedmachine thickness along outlet perimeter pixels is 1.047*hydrostatic (see thickness_inversion_uncertainty.m)    
th_bm(mask_bm==0) = th_bm_inv(mask_bm==0);

% Overwrite BedMachine thickness where AERODEM data are available: 
aerodem = load('aerodem_clean.mat'); 
was_floating = aerodem.sfzk*9.34 < th_bm_inv; 
th_bm(aerodem.marine_retreat & was_floating) = 9.34*aerodem.sfzk(aerodem.marine_retreat & was_floating);
th_bm(aerodem.marine_retreat & ~was_floating) = aerodem.sfzk(aerodem.marine_retreat & ~was_floating) - bed_bm(aerodem.marine_retreat & ~was_floating);
errbed_bm(aerodem.marine_retreat & was_floating) = 93.4; % 10 m surface accuracy error, inverted for thickness.   
errbed_bm(aerodem.marine_retreat & ~was_floating) = hypot(10,errbed_bm(aerodem.marine_retreat & ~was_floating));

% Interpolate from BedMachine grid to our grid: 
th = interp2(x_bm,y_bm,th_bm,X,Y); 
th_error = interp2(x_bm,y_bm,errbed_bm,X,Y); 
maskbm = interp2(x_bm,y_bm,mask_bm,X,Y,'nearest'); 
maskbm(maskbm==4) = 0; 
maskbm(isnan(maskbm)) = 0; 

fn = 'GimpOceanMask_90m_2015_v1.2.tif'; 
[mx,my]=pixcenters(geotiffinfo(fn));
ocean = interp2(mx,my,logical(imread(fn)),X,Y,'nearest',1); 
ocean = ~imfill(~ocean,8,'holes'); 

fn = 'GimpIceMask_90m_2015_v1.2.tif'; 
ice = interp2(mx,my,logical(imread(fn)),X,Y,'nearest',0); 

% Override masks where BedMachine says it's rock, then combine to create a rock mask: 
ocean(maskbm==1) = false; 
ice(maskbm==1) = false;
rock = ~ocean & ~ice; 

th(rock) = 0; 
thickness_source = 3*ones(size(th),'uint8');          % bed inversion
thickness_source(th>0 & ismember(maskbm,[2 3])) = 1;  % bedmachine
thickness_source(interp2(x_bm,y_bm,aerodem.marine_retreat,X,Y,'nearest',0)==1) = 2; % aerodem 
thickness_source(rock) = 0; 
th_error(rock) = 0; 

S = shaperead('/Users/cgreene/Documents/data/basins/doi_10.7280_D1WT11__v1/Greenland_Basins_PS_v1.4.2.shp'); 

clear th_bm x_bm y_bm bed_bm th_bm_inv mask_bm mx my maskbm fn errbed_bm aerodem was_floating
disp 'data loaded'

%% Glacier catchment masking 
% Catchments 226-230 take a while due to some polygons with a lot of vertices. 

try 
   % If I've already done this, there's no need to do it again, so just load the old data. 
   catchment = permute(ncread('/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_extruded_velocity_and_thickness_2022-12-08.nc','catchment'),[2 1]); 

catch
   catchment = zeros(length(y),length(x),'int16'); 
   
   for k = 1:length(S)
      catchment(inpolygon_map(double(x),double(y),S(k).X,S(k).Y)) = k; 
      disp(['Catchment ',num2str(k),' of 260'])
   end
end

disp 'initial catchment done'

%% Standardize gridded datasets 

vx(rock) = 0; 
vy(rock) = 0; 

%%

% Holes are in the interior of the ice sheet, and are surrounded by valid
% ITS_LIVE measurements. Fill them with regionfill and not model data. 
holes = isnanhole(vx); 

v_source = 4*ones(size(vx),'uint8'); % 4 = extrapolated willy nilly 
v_source(isfinite(vx)) = 1;          % 1 = observations
v_source(rock) = 0;                  % 0 = rock. 
v_source(holes) = 2;                 % 2 = interpolated from observations 
                                     % 3 = extrapolated along paleo ice sheet model flowlines (below) 
                                     
vx = regionfill(vx,holes);
vy = regionfill(vy,holes);

clear holes 

%% Modeled flow velocities 
% Paleo ice sheet data from Josh Cuzzone 

load('/Users/cgreene/Documents/data/velocity/GrIS_Paleo_IceFlow_For_Chad_v3/lat_GrIS_JPL.mat')
load('/Users/cgreene/Documents/data/velocity/GrIS_Paleo_IceFlow_For_Chad_v3/lon_GrIS_JPL.mat')
load('/Users/cgreene/Documents/data/velocity/GrIS_Paleo_IceFlow_For_Chad_v3/Vx_GrIS_JPL.mat')
load('/Users/cgreene/Documents/data/velocity/GrIS_Paleo_IceFlow_For_Chad_v3/Vy_GrIS_JPL.mat')

v_vx = flipud(permute(v_vx,[2 3 1])); 
v_vy = flipud(permute(v_vy,[2 3 1])); 
lat = flipud(lat); 

% Combine model velocities. 
% The first timestep is paleo, the last one is present day. So work
% backwards, starting with the ice sheet that's most like today, and fill
% in the missing pieces with older and older ice as needed: 
vx_tmp = nan(size(v_vx,1),size(v_vx,2)); 
vy_tmp = vx_tmp; 
for k = size(v_vx,3):-1:1
   tmp = v_vx(:,:,k);
   vx_tmp(isnan(vx_tmp)) = tmp(isnan(vx_tmp)); 
   tmp = v_vy(:,:,k);
   vy_tmp(isnan(vy_tmp)) = tmp(isnan(vy_tmp)); 
end

% Get model data on ITS_LIVE grid: 
[Lat,Lon] = psn2ll(X,Y);
vx_mod = interp2(lon,lat,vx_tmp,Lon,Lat);
vy_mod = interp2(lon,lat,vy_tmp,Lon,Lat);

% Fill itslive velocities with model data: 
v_source(isnan(vx) & isfinite(vx_mod)) = 3; 

vx(v_source==3) = vx_mod(v_source==3); 
vy(v_source==3) = vy_mod(v_source==3); 

% If infilled holes are actually just tight rock-bound fjords that need modeled velocities: 
fix = v_source==2 & hypot(vx,vy)<hypot(vx_mod,vy_mod); 
vx(fix) = vx_mod(fix); 
vy(fix) = vy_mod(fix); 
v_source(fix) = 3; 

clear v_vx v_vy lon lat Lon Lat vx_mod vy_mod fix 
disp 'modeled velocities complete'

%% Extrude observed vel through modeled vel region 

% Get the boundary of where we have observations: 
has = bwselect(imfill(v_source<=2,'holes'),floor(size(rock,2)/2),floor(size(rock,1)/2)); 

v = hypot(vx,vy); 

% Use the perimeter of the velocity measurements as the seed locations: 
outlets = bwperim(has) & ~rock;
outlets(bwperim(true(size(rock)))) = false; 
outlets(v==0) = false; 
[row,col] = find(outlets); 
XY = stream2(double(x),double(y),vx,vy,double(X(outlets)),double(Y(outlets)),[.2 10000]);

% Loop through each seed location to extrapolate terminus thickness along flowlines:  
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = v(row(k),col(k)); 
end

% Concatenate the cells: 
M = cell2mat(XY(:)); % 
V = cell2mat(V(:)); 

% Grid up the streamline data:  
isf = isfinite(V(:,1)) & isfinite(M(:,1)); 
vg = gridbin(M(isf,1),M(isf,2),V(isf,1),x,y); %

% Fill gridded obs with directly observed obs or interpolations: 
vg(v_source<=2) = v(v_source<=2);

% Only fill in the holes between streamlines
fll = isnanhole(vg);

% But first fill in the nans that aren't holes: 
vg(isnan(vg) & ~fll) = v(isnan(vg) & ~fll); 
vg(rock) = 0; 

% Now just fill the streamline holes: 
vg = regionfill(vg,fll); 

% To prevent near-zero velocity streaks, replace low streaky velocities with modeled (generally conservative) velocities:   
vg(v>vg & v_source==3) = v(v>vg & v_source==3); 

% Overwrite vx and vy as original unit vectors scaled by vg:  
vx = vg.*vx./v; 
vy = vg.*vy./v; 
vx(rock) = 0; 
vy(rock) = 0; 

clear vg fll v has isf outlets 
disp 'first velocity extrapolation complete'

%% Extrapolate remaining Velocity 

skip = 4; 
xr = x(1:skip:end); 
yr = y(1:skip:end); 

% Downsample by a factor of 4 (to ~500 m resolution) to make inpaint_nans tractable):  
Vxr = inpaint_nans(th(1:skip:end,1:skip:end).*vx(1:skip:end,1:skip:end),4); % takes 1 minute with skip=4 
Vyr = inpaint_nans(th(1:skip:end,1:skip:end).*vy(1:skip:end,1:skip:end),4); % takes 1 minute with skip=4 
Vr = hypot(Vxr,Vyr); 

v = hypot(vx,vy); 

% Use the perimeter of the velocity measurements as the seed locations: 
outlets = bwperim((isfinite(v) | rock)) & ~rock;
outlets(bwperim(true(size(rock)))) = false; 
[row,col] = find(outlets); 
XY = stream2(double(xr),double(yr),Vxr,Vyr,double(X(outlets)),double(Y(outlets)),[0.2 25e3]);

% Loop through each seed location to extrapolate terminus velocity and ice shelf name along flowlines:  
V = XY; 
for k = 1:length(V) 
   V{k}(:,1) = v(row(k),col(k)); 
end

% Concatenate the cells: 
M = cell2mat(XY(:)); % 
V = cell2mat(V(:)); 

% Grid up the streamline data:  
isf = isfinite(V(:,1)) & isfinite(M(:,1)); 
vg = gridbin(M(isf,1),M(isf,2),V(isf,1),x,y); %

% Wherever gridbin gives us data, overwrite tmpvx, tmpvy:  
isf = isfinite(vg) & isnan(v) & ~rock;
vx(isf) = interp2(xr,yr,Vxr./Vr,X(isf),Y(isf)).*vg(isf); 
vy(isf) = interp2(xr,yr,Vyr./Vr,X(isf),Y(isf)).*vg(isf); 

% Make sure rock and non-greenland land are zeros: 
vx(rock) = 0; % rock or non-Greenland land
vy(rock) = 0; 

% Fill holes: 
holes = isnanhole(vx); 
vx = regionfill(vx,holes); 
vy = regionfill(vy,holes); 

% Fill anything that's left over with zero values: 
vx(isnan(vx)) = 0; 
vy(isnan(vy)) = 0; 

clear Vxr Vyr XY Vr xr yr 

disp 'second velocity extrapolation complete'

%% Catchment extrapolation 

rockholes = imfill(rock,8,'holes') & ~rock;

% Use the perimeter of the velocity measurements as the seed locations: 
outlets = bwperim(catchment>0 | rock) & ~rock;
outlets(bwperim(true(size(rock)))) = false; 
outlets(v==0) = false; 
[row,col] = find(outlets); 

for k = 1:260
   % Streamlines from the outlet pixels of this glacier's catchment: 
   XY = stream2(double(x),double(y),vx,vy,double(X(outlets & catchment==k)),double(Y(outlets & catchment==k)),[.2 10000]);

   M = cell2mat(XY(:)); % 
   if ~isempty(M)
      isf = isfinite(M(:,1)); 
   
      % Grid up the down-streamlines: 
      tmp = gridbin(M(isf,1),M(isf,2),true(size(M(isf,1))),x,y,@any); 
   
      % Overwrite the catchment numbers where the original dataset says it's ocean   
      catchment(imfill(tmp | catchment==k | rock,4,'holes') & catchment==0 & ~rock & ~rockholes) = k; 
   end

   disp(['Catchment extrapolation ',num2str(k),' of 260.'])
end

% Dilate each basin by 5 pixel radius (600 m) 
st = strel('disk',5,0); 

for kk = 1:10
   for k = 1:260
      catchment(catchment==0 & imdilate(catchment==k & ~rock,st) & ~rock & ~rockholes) = k; 
   end
   kk

end

catchment(catchment==0) = 261; % "other"
catchment(rock) = 0; 

figure
imagescn(x,y,catchment)
axis image off
bedmachine('gl','k','greenland')
rng('default')
colormap(rand(1e3,3))

%% Save data 

filename = ['/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_extruded_velocity_and_thickness_',datestr(now,'yyyy-mm-dd'),'.nc']; 

disp 'saving data now...'
th(isnan(th)) = 32767;
thickness = uint16(ceil(th)); % NC_USHORT
thickness_error = uint16(ceil(th_error)); 
vx = int16(round(vx)); % NC_SHORT
vy = int16(round(vy)); % NC_SHORT

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(filename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Greenland extruded velocity and thickness');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','This dataset imagines a Greenland Ice Sheet without calving fronts. Velocity and thickness observations are extrapolated beyond present-day extents, along paleo flow directions where available (and totally made up beyond that), to fill the entire domain. Thickness and velocity extrapolations are constant, and therefore mass rates are somewhat conserved, except that flow may converge or diverge, in which case volume is not conserved. Close to the present-day calving fronts, however, this dataset provides a reasonable first-order approximation of thickness and velocity beyond their observed extents.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date_created',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Institution','NASA Jet Propulsion Laboratory (JPL), California Institute of Technology');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Please cite this dataset! Check the GitHub page for citation info.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'GitHub','https://github.com/chadagreene/greenland-coastlines');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'MATLAB_script','This dataset was generated by greenland_extrude_velocity_and_thickness_v2.m.');

% 2. Define dimensions
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'false_easting',0);
netcdf.putAtt(ncid,mapping_var_id,'false_northing',0);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid,mapping_var_id,'inverse_flattening',298.257224);
netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',90);
netcdf.putAtt(ncid,mapping_var_id,'scale_factor_at_projection_origin',1);
netcdf.putAtt(ncid,mapping_var_id,'semi_major_axis',6378137);
netcdf.putAtt(ncid,mapping_var_id,'spatial_epsg',3413);
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj','+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs');
netcdf.putAtt(ncid,mapping_var_id,'spatial_ref','PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",SOUTH],AXIS["Northing",SOUTH],AUTHORITY["EPSG","3413"]]');
netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',70);
netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',-45);

% Define x 
x_id     = netcdf.defDim(ncid,'x',length(x));
x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',x_id);
netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate');
netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,x_var_id,'units',        'meter');

% Define y
y_id     = netcdf.defDim(ncid,'y',length(y));
y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',y_id);
netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate');
netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,y_var_id,'units',        'meter');

% Define vx
vx_var_id = netcdf.defVar(ncid,'vx','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid, vx_var_id,'long_name','projected x component of velocity');
netcdf.putAtt(ncid, vx_var_id,'standard_name','vx');
netcdf.putAtt(ncid, vx_var_id,'units',    'm/yr');
netcdf.putAtt(ncid, vx_var_id,'grid_mapping', 'mapping');

% Define vy
vy_var_id = netcdf.defVar(ncid,'vy','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid, vy_var_id,'long_name','projected y component of velocity');
netcdf.putAtt(ncid, vy_var_id,'standard_name','vy');
netcdf.putAtt(ncid, vy_var_id,'units',    'm/yr');
netcdf.putAtt(ncid, vy_var_id,'grid_mapping', 'mapping');

% Define v_source
v_source_var_id = netcdf.defVar(ncid,'v_source','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid, v_source_var_id,'long_name',    '0=rock; 1=observations; 2=interpolated observations; 3=extrapolated observations following paleo model flowlines; 4=total extrapolation');
netcdf.putAtt(ncid, v_source_var_id,'grid_mapping', 'mapping');

% Define thickness
thick_var_id = netcdf.defVar(ncid,'thickness','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid,thick_var_id,'long_name','ice thickness');
netcdf.putAtt(ncid,thick_var_id,'standard_name','ice_thickness');
netcdf.putAtt(ncid,thick_var_id,'_FillValue',int16(Inf));
netcdf.putAtt(ncid,thick_var_id,'units',    'meters');
netcdf.putAtt(ncid,thick_var_id,'grid_mapping', 'mapping');

% Define thickness_error
thick_err_var_id = netcdf.defVar(ncid,'thickness_error','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid,thick_err_var_id,'long_name','ice thickness error taken as errbed for present-day grounded ice and hypot(1.12*1.047*errbed,32.7) aka the root-sum-square of bed elevation error and the error associated with the inversion.');
netcdf.putAtt(ncid,thick_err_var_id,'standard_name','ice_thickness_error');
netcdf.putAtt(ncid,thick_err_var_id,'_FillValue',int16(Inf));
netcdf.putAtt(ncid,thick_err_var_id,'units',    'meters');
netcdf.putAtt(ncid,thick_err_var_id,'grid_mapping', 'mapping');

% Define thickness_source
thickness_source_var_id = netcdf.defVar(ncid,'thickness_source','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid, thickness_source_var_id,'long_name',    '0=rock; 1=BedMachine v5; 2=AERODEM; 3=bed inversion');
netcdf.putAtt(ncid, thickness_source_var_id,'grid_mapping', 'mapping');

% Define catchment
catchment_var_id = netcdf.defVar(ncid,'catchment','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid, catchment_var_id,'long_name',    'Glacier catchments extrapolated from Mouginot, Jeremie; Rignot, Eric (2019), Glacier catchments/basins for the Greenland Ice Sheet, Dryad, Dataset, https://doi.org/10.7280/D1WT11.');
netcdf.putAtt(ncid, catchment_var_id,'name_file',    'Glacier catchments names are provided in glacier_catchment_names_Mouginot.xlsx.');
netcdf.putAtt(ncid, catchment_var_id,'grid_mapping', 'mapping');
netcdf.putAtt(ncid,catchment_var_id,'_FillValue',int16(Inf));

% Compress and stop variable definition
netcdf.defVarDeflate(ncid,vx_var_id,true,true,9);
netcdf.defVarDeflate(ncid,vy_var_id,true,true,9);
netcdf.defVarDeflate(ncid,v_source_var_id,true,true,9);
netcdf.defVarDeflate(ncid,thick_var_id,true,true,9);
netcdf.defVarDeflate(ncid,thick_err_var_id,true,true,9);
netcdf.defVarDeflate(ncid,thickness_source_var_id,true,true,9);
netcdf.defVarDeflate(ncid,catchment_var_id,true,true,9);
netcdf.endDef(ncid);

%3. Place data
netcdf.putVar(ncid,x_var_id,x);
netcdf.putVar(ncid,y_var_id,y);
netcdf.putVar(ncid,vx_var_id,ipermute(vx,[2 1]));
netcdf.putVar(ncid,vy_var_id,ipermute(vy,[2 1]));
netcdf.putVar(ncid,v_source_var_id,ipermute(v_source,[2 1]));
netcdf.putVar(ncid,thick_var_id,ipermute(thickness,[2 1]));
netcdf.putVar(ncid,thick_err_var_id,ipermute(thickness_error,[2 1]));
netcdf.putVar(ncid,thickness_source_var_id,ipermute(thickness_source,[2 1]));
netcdf.putVar(ncid,catchment_var_id,ipermute(catchment,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp 'all done'

%% Subfunctions 
function holes = isnanhole(bw)
   holes = imfill(isfinite(bw),8,'holes') & ~isfinite(bw);
end




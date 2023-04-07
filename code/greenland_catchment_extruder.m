% After running greenland_extrude_velocity_and_thickness_v2.m, I discovered
% some problems with the catchment definition and extrusion. This script 
% fixes that, and rewrites the extruded thickness, velocity, and catchment
% file. 
% 
% This script takes several hours to run. 
% 
% Chad Greene NASA/JPL, April 2023. 

fn = 'greenland_extruded_velocity_and_thickness_2022-12-15.nc'; 
x = ncread(fn,'x'); 
y = ncread(fn,'y'); 
vx = permute(ncread(fn,'vx'),[2 1]); 
vy = permute(ncread(fn,'vy'),[2 1]); 
rock = logical(permute(ncread('greenland_ice_masks_1972-2022.nc','rock'),[2 1])); 


[X,Y] = meshgrid(x,y); 
ice0 = bedmachine_interp('mask',X,Y,'greenland')==2; 

S = shaperead('/Users/cgreene/Documents/data/basins/doi_10.7280_D1WT11__v1/Greenland_Basins_PS_v1.4.2.shp'); 


catchment = zeros(length(y),length(x),'int16'); 

for k = 1:length(S)
    catchment(inpolygon_map(double(x),double(y),S(k).X,S(k).Y)) = k; 
    disp(['Catchment ',num2str(k),' of 260'])
end

catchment0 = catchment; 

%% Extrude 

outlets = imdilate(~(catchment>0 | rock ),true(9)); 
outlets(imfill(rock,'holes')) = false; 
outlets(vx==0 & vy==0) = false; 
outlets(catchment==0) = false; 
outlets(bwdist(~imfill(catchment0>0,'holes'))>15) = false; 

% figure
% imagescn(x,y,outlets)
% bedmachine('gl','greenland') 

%% 

%rockholes = imfill(rock,8,'holes') & ~rock;

% % Use the perimeter of the velocity measurements as the seed locations: 
% outlets = bwperim(catchment>0 | rock) & ~rock;
% outlets(bwperim(true(size(rock)))) = false; 
% outlets(v==0) = false; 
[row,col] = find(outlets); 

for k = 1:260
   % Streamlines from the outlet pixels of this glacier's catchment: 
   XY = stream2(double(x),double(y),vx,vy,double(X(outlets & catchment==k)),double(Y(outlets & catchment==k)),[.2 30000]);

   M = cell2mat(XY(:)); % 
   if ~isempty(M)
      isf = isfinite(M(:,1)); 
   
      % Grid up the down-streamlines: 
      tmp = gridbin(M(isf,1),M(isf,2),true(size(M(isf,1))),x,y,@any); 
   
      % Overwrite the catchment numbers where the original dataset says it's ocean   
      %catchment(imfill(tmp | catchment==k | rock,4,'holes') & catchment==0 & ~rock & ~rockholes) = k; 
      catchment(tmp & ~rock & catchment==0) = k; 
   end

   disp(['Catchment extrapolation ',num2str(k),' of 260.'])
end

clear vx vy X Y
% save test.mat

%% 
disp 'filling holes '
for k = 1:260 
    catchment(imfill(catchment==k,'holes') & ~rock & catchment==0) = k; 
    k
end   

%% 

disp dilating 
rockfill = imfill(rock,8,'holes'); 

st = strel('disk',9);

for kn = 1:10
    for k = 1:260 

        tmp = imfill(imdilate(catchment==k,st),'holes'); 
    
        catchment(tmp & ~rockfill & catchment==0) = k; 
    end                 
    kn
end

%%

catchment(~rockfill & inpolygon_map(double(x),double(y),[509452      509162      523066      534556      532529      522583],[ -1061541    -1080562    -1092922    -1082107    -1058548    -1049375])) = 59; 


catchment(catchment==0) = 261; % "other"
catchment(rock) = 0; 

%%
T = load('terminus_data_densified_2023-01-09.mat');
spacing_m = 120/5; % along-path spacing of the densified terminus data is 5 points per grid cell. 
t = datenum(datetime(1972,9,15):calmonths(1):datetime(2022,2,15));
[yr,mo,~] = datevec(t); 

[yr_data,mo_data,~] = datevec(double(T.t)); 

catchmenti = interp2(x,y,catchment,T.x,T.y,'nearest'); 

obs_km = zeros(numel(t),261); 
disp 'observations' 
for k=1:length(t)
%     for catchk = 1:261
%         obs_km(k,catchk) = sum(catchmenti==catchk & yr_data==yr(k) & mo_data==mo(k))*spacing_m/1000; 
%     end
    ind = yr_data==yr(k) & mo_data==mo(k);
    obs_km(k,:) = histcounts(catchmenti(ind),.5:261.5)*spacing_m/1000;
    k
end

% save('greenland_catchments_2023-04-06.mat','x','y','t','obs_km','names','catchment')
%%

oldfilename = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_extruded_velocity_and_thickness_2022-12-15.nc';

x = ncread(oldfilename,'x'); 
y = ncread(oldfilename,'y'); 
vx = ncread(oldfilename,'vx'); 
vy = ncread(oldfilename,'vy'); 
v_source = ncread(oldfilename,'v_source'); 
thickness = ncread(oldfilename,'thickness'); 
thickness_error = ncread(oldfilename,'thickness_error'); 
thickness_source = ncread(oldfilename,'thickness_source'); 


thickness(isnan(thickness)) = 32767;
thickness = uint16(ceil(thickness)); % NC_USHORT
thickness_error = uint16(ceil(thickness_error)); 

%%

newfilename = ['/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_extruded_velocity_and_thickness_',datestr(now,'yyyy-mm-dd'),'.nc']; 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Greenland extruded velocity and thickness');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','This dataset imagines a Greenland Ice Sheet without calving fronts. Velocity and thickness observations are extrapolated beyond present-day extents, along paleo flow directions where available (and totally made up beyond that), to fill the entire domain. Thickness and velocity extrapolations are constant, and therefore mass rates are somewhat conserved, except that flow may converge or diverge, in which case volume is not conserved. Close to the present-day calving fronts, however, this dataset provides a reasonable first-order approximation of thickness and velocity beyond their observed extents.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date_created',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Institution','NASA Jet Propulsion Laboratory (JPL), California Institute of Technology');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Please cite this dataset! Check the GitHub page for citation info.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'GitHub','https://github.com/chadagreene/greenland-icemask');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'MATLAB_script','This dataset was generated by greenland_extrude_velocity_and_thickness_v2.m, then edited by greenland_catchment_extruder.m.');

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
%netcdf.putAtt(ncid,catchment_var_id,'_FillValue',int16(Inf));

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
netcdf.putVar(ncid,vx_var_id,vx);
netcdf.putVar(ncid,vy_var_id,vy);
netcdf.putVar(ncid,v_source_var_id,v_source);
netcdf.putVar(ncid,thick_var_id,thickness);
netcdf.putVar(ncid,thick_err_var_id,thickness_error);
netcdf.putVar(ncid,thickness_source_var_id,thickness_source);
netcdf.putVar(ncid,catchment_var_id,ipermute(catchment,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp 'all done'



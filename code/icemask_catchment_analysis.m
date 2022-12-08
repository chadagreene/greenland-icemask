

% This script tallies up mass and area changes per glacier basin 
% using Greenland ice masks. 
% 
% The script is super slow, but who cares, just give it like 14 hours to
% run. 
% 
% Chad A. Greene, NASA/JPL, November 2022. 

%% Load data 

extruded_filename = 'greenland_extruded_velocity_and_thickness_2022-11-29.nc'; 
icemask_filename = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_monthly_ice_masks_2022-11-18_clean.nc';

x = double(ncread(extruded_filename,'x')); 
y = double(ncread(extruded_filename,'y')); 

rock = permute(ncread(extruded_filename,'v_source'),[2 1])==0; 
th = double(permute(ncread(extruded_filename,'thickness'),[2 1])); 
th_err = double(permute(ncread(extruded_filename,'thickness_error'),[2 1])); 
catchment = double(permute(ncread(extruded_filename,'catchment'),[2 1])); 

load('icemask_catchment_analysis.mat','names') % it's recursive, I know, but these are the Mouginot catchment names. 

%%

[X,Y] = meshgrid(x,y); 
[Lat,~] = psn2ll(X,Y); 
Area = (diff(x(1:2)) ./ psndistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  

Mass = Area.*th.*917*1e-12; 
Mass_hi = Area.*(th+th_err).*917*1e-12; 
clear Lat X Y th

%%

readme = 'Ice time series for Mouginot catchments 1-260 plus "other". Area in m^2 and Mass in Gt. Created by icemask_catchment_analysis.m.'; 

t = ncdateread(icemask_filename,'time');

A_ts = nan(length(t),max(catchment(:))); 
M_ts = A_ts; 
M_hi_ts = M_ts; 

for k = length(t):-1:1
    ice_tmp = permute(logical(ncread(icemask_filename,'ice',[1 1 k],[Inf Inf 1])),[2 1]);
    
    catchment2 = cube2rect(catchment,ice_tmp); 
    Area2 = cube2rect(Area,ice_tmp); 
    Mass2 = cube2rect(Mass,ice_tmp); 
    Mass_hi2 = cube2rect(Mass_hi,ice_tmp); 
    for kb = 1:261
        A_ts(k,kb) = sum(Area2(catchment2==kb)); 
        M_ts(k,kb) = sum(Mass2(catchment2==kb)); 
        M_hi_ts(k,kb) = sum(Mass_hi2(catchment2==kb)); 
    end
    disp([datestr(now),' month ',num2str(k),' of ',num2str(length(t)),' complete.'])

    %save('/Users/cgreene/Documents/GitHub/greenland-coastlines/data/icemask_catchment_analysis_clean.mat','t','A_ts','M_ts','M_hi_ts','readme','names','icemask_filename','extruded_filename')
end

%%


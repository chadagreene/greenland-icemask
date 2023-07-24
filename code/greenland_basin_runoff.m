
load('ice_catchment_analysis_area_mass_totals.mat','names')
load('greenland_ice_masks_1972-2022_v1_ice_sum.mat','ice_sum')
x = ncread('greenland_ice_masks_1972-2022_v1.nc','x'); 
y = ncread('greenland_ice_masks_1972-2022_v1.nc','y'); 

catchment = permute(ncread('greenland_ice_masks_1972-2022_v1.nc','catchment'),[2 1]); 


[X,Y] = meshgrid(x,y); 

[Lat,~] = psn2ll(X,Y); 
A = (diff(x(1:2)) ./ psndistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  

fn = '/Users/cgreene/Documents/data/SMB/gsfc_fdm_smb_v1_2_1_gris_Dec22.nc'; 
x_gsfc = ncread(fn,'x'); 
y_gsfc = ncread(fn,'y'); 
t_gsfc = ncread(fn,'time'); 
runoff = ncread(fn,'Ru');

grid_cell_area = 12500^2; % m^2
rho_ice = 917;            
kg_to_Gt = 1e-12; % multiplier for converting kg to Gt 
s1 = squeeze(sum(sum(runoff,'omitnan'),'omitnan'))*grid_cell_area*rho_ice*1e-12;

t_range = t_gsfc>=2000 & t_gsfc<2021; 
mean_runoff_m_ice_per_year = mean(runoff,3,'omitnan'); % grid of mean rate of runoff.

total_runoff_Gt_per_year = sum(mean_runoff_m_ice_per_year(:),'omitnan') * grid_cell_area * rho_ice * kg_to_Gt

t_range = t_gsfc>=2000 & t_gsfc<2021; 

runoff_annual = interp2(single(x_gsfc),single(y_gsfc),mean(runoff(:,:,t_range),3),X,Y) .*A * rho_ice * kg_to_Gt;

clear runoff t_gsfc x_gsfc y_gsfc fn

runoff_basin = nan(size(names)); 
for k = 1:261
    runoff_basin(k)  = sum(runoff_annual(catchment==k & ice_sum==594)); 
end

basin_area = nan(size(names)); 
for k = 1:261
    basin_area(k) = sum(A(ice_sum==594 & catchment==k)); 
end

runoff_basin_mmwe = runoff_basin * 1e12 ./basin_area;

save('/Users/cgreene/Documents/GitHub/greenland-icemask/data/greenland_basin_runoff.mat','runoff_basin','basin_area','runoff_basin_mmwe','names')

% Use TopoToolbox to calculate the sill depth for each basin. 
% Chad Greene July 2023. 


fn = 'greenland_extruded_velocity_and_thickness_2023-04-06.nc';
load('greenland_ice_masks_1972-2022_v1_ice_sum.mat')

th = permute(ncread(fn,'thickness'),[2 1]); 
catchment = permute(ncread(fn,'catchment'),[2 1]); 

[X,Y] = meshgrid(x,y); 
bed = bedmachine_interp('bed',X,Y,'greenland'); 

%% 

Z = bed + th;
Z(ice_sum<594) = bed(ice_sum<594);
Zg  = GRIDobj(x,y,Z); % topotoolbox 
Zgf  = fillsinks(Zg);

bed_filled = Zgf.Z; 

sill_depth  = nan(261,1); 
for k = 1:261
    msk = catchment==k & ice_sum<594 & ice_sum>0; 
    if any(msk,'all')
        sill_depth(k) = min(bed_filled(msk)); 
    end
end

%% 

readme = 'sill depth corresponding to each of Mouginots greenland basins. Saved by sill_depth_calculator.m.'; 

%save('/Users/cgreene/Documents/GitHub/greenland-icemask/data/sill_depth.mat','sill_depth','readme') 

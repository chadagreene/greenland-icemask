

load('ice_catchment_analysis_area_mass_totals.mat','names')
load('greenland_ice_masks_1972-2022_v1_ice_sum.mat','ice_sum')
rock = permute(logical(ncread('greenland_ice_masks_1972-2022_v1.nc','rock')),[2 1]); 

fn = 'greenland_extruded_velocity_and_thickness_2023-04-06.nc';
vx = double(permute(ncread(fn,'vx'),[2 1])); 
vy = double(permute(ncread(fn,'vy'),[2 1])); 
th = permute(ncread(fn,'thickness'),[2 1]); 
catchment = permute(ncread(fn,'catchment'),[2 1]); 
x = double(ncread(fn,'x')); 
y = double(ncread(fn,'y')); 

[X,Y] = meshgrid(x,y); 


%%

outlets = bwperim(ice_sum==594 | rock) & ~rock & hypot(vx,vy)>0 & th>0; 

step = 0.2; % Number of streamline steps per grid cell.  
reach = 8e3; % meters of reach in each direction from measured termini
Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 
       

bed_slope = nan(261,1); 
sfz_slope = bed_slope; 

for catchk = 1:261
    xo = X(outlets & catchment==catchk); 
    yo = Y(outlets & catchment==catchk); 
    tho = th(outlets & catchment==catchk); 
    ds = stream2(x,y,-vx,-vy,xo,yo,[step Npts]);
    
    bed_slopekk = nan(length(ds),1); 
    sfz_slopekk = bed_slopekk;
    for kk = 1:length(ds)
        dkk = pathdistpsn(ds{kk}(:,1),ds{kk}(:,2)); 
        ind = dkk<=5000; 
        dkk = dkk(ind); 
        if dkk(end)>4800
            bedkk = bedmachine_interp('bed',ds{kk}(ind,1),ds{kk}(ind,2),'greenland'); 
    
            sfzkk = bedkk + interp2(x,y,th,ds{kk}(ind,1),ds{kk}(ind,2));
    
            pv = polyfit(-dkk,bedkk,1);
            bed_slopekk(kk) = pv(1); 
    
            pv = polyfit(-dkk,sfzkk,1);
            sfz_slopekk(kk) = pv(1); 
        end
    end
       
    isf = isfinite(bed_slopekk);
    if any(isf)
        bed_slope(catchk) = atand(wmean(bed_slopekk(isf),tho(isf)));
    end

    isf = isfinite(sfz_slopekk);
    if any(isf)
        sfz_slope(catchk) = atand(wmean(sfz_slopekk(isf),tho(isf)));
    end

    catchk
end

%% 

readme = 'thickness-weighted mean bed and surface slopes (calculated by least squares) from all perimeter grid cells in each catchment, within 5 km of the most retreated position of each glacier. Calculated by catchment_bed_slopes.m';

%save('/Users/cgreene/Documents/GitHub/greenland-icemask/data/catchment_bed_slopes.mat','sfz_slope','bed_slope','names','readme')


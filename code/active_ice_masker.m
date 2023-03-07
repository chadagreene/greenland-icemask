

icemask_filename = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_monthly_ice_masks_2023-02-22_compressed.nc'; 

t = ncdateread(icemask_filename,'time');
x = ncread(icemask_filename,'x');
y = ncread(icemask_filename,'y');


t1 = find(t>=datetime('june 15, 1985'),1,'first'); 
t2 = find(t>=datetime('june 15, 2021'),1,'first'); 

always_ice = logical(ncread(icemask_filename,'ice',[1 1 t1],[Inf Inf 1]));
never_ice = ~always_ice; 

for k = (t1+1):t2
    ice_tmp = logical(ncread(icemask_filename,'ice',[1 1 k],[Inf Inf 1]));
    
    always_ice(~ice_tmp) = false; 
    never_ice(ice_tmp) = false; 
    k
end

always_ice = permute(always_ice,[2 1]); 
never_ice = permute(never_ice,[2 1]); 
active_ice = ~always_ice & ~never_ice;

readme = 'a quick mask of pixels that changed between ice and ocean anytime between June 1985 and June 2021.';

%save('/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/active_ice_mask_2023-01-24.mat','x','y','active_ice','always_ice','never_ice','readme')

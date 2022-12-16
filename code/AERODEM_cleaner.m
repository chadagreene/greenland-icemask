% After creating a mosaic from Korsgaard's AERODEM, run this to clean 
% up the silly data bits. The whole script takes 10 seconds to run. 
% 
% Chad A. Greene, NASA/JPL, December 2022.


%% Load data 

filename = 'AERODEM_Korsgaard_mosaic.nc'; 

x = ncread(filename,'x'); 
y = ncread(filename,'y'); 
sfzk = permute(ncread(filename,'surface'),[2 1]);
reliab = permute(ncread(filename,'reliability_mask'),[2 1]);

sfzb = bedmachine_data('surface','greenland'); 
rock = bedmachine_data('mask','greenland')==1; 
bed = bedmachine_data('bed','greenland'); 

%% Adjust data to taste

sfzk(reliab<40) = nan; % expected low reliability 
sfzk(sfzk<10) = nan; % Gets rid of sea ice and other inconsequential bits 

holes = isnanhole(sfzk); % isnanhole is a function at the bottom of this script. 

sfzk = regionfill(sfzk-sfzb,holes) + sfzb; % fills in the holes  

sfzk(~bwareafilt(isfinite(sfzk),[1e3 Inf])) = nan; % eliminates small icebergs and what-have-you.  

figure
imagescn(x,y,sfzk-sfzb)
axis image off
caxis([-1 1]*100)
cmocean bal

%%

% Marine retreat has occurred where bed is below sea level, and 1980's surface is above sea level and BedMachine says there's no ice at all.  
marine_retreat = sfzk>0 & sfzb==0 & bed<0; 

figure
imagescn(x,y,marine_retreat)
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')

%%

readme = 'Surface elevations, wrt geoid, from AERODEM, with bad data removed and missing holes filled in. Marine retreat has occurred where bed is below sea level, and 1980s surface is above sea level and BedMachine says there is no ice at all. ';

% save('/Users/cgreene/Documents/data/DEMs/korsgaard_g150_aerodem/aerodem_clean.mat','sfzk','marine_retreat','x','y','readme')

disp 'done' 
%% Subfunctions 
function holes = isnanhole(bw)
   holes = imfill(isfinite(bw),4,'holes') & ~isfinite(bw);
end
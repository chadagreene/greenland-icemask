% This script compiles Measures Greenland annual terminus positions v2 
% (Joughin et al., 2021) https://nsidc.org/data/NSIDC-0642) into a single .mat file. 
% The reason for this script is because the data when downloaded are spread
% across more than 400 files (shapefile and all of the ancillary files that
% come with shapefiles), and each year of data are in arbitrarily named
% folders. 
% 
% Before running this script, manually move all shapefiles and their
% ancillary files into one folder. 
% 
% Chad Greene, July 2022. 

f = dir('/Users/cgreene/Documents/data/coastlines/measures_greenland_terminus_v2/termini*_v02.0.shp');

S = shaperead(['/Users/cgreene/Documents/data/coastlines/measures_greenland_terminus_v2/',f(1).name]);

for k = 2:length(f) 
   S = cat(1,S,shaperead(['/Users/cgreene/Documents/data/coastlines/measures_greenland_terminus_v2/',f(k).name]));
end
   
readme = 'This is Measures data (Joughin et al., 2021 https://nsidc.org/data/NSIDC-0642) repackaged from 49 shapefiles to a single .mat file by measures_terminuspicks_compiler.m.';

%save('/Users/cgreene/Documents/GitHub/greenland-coastlines/data/measures_greenland_terminus_picks_v2.mat','S','readme')



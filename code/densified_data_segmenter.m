% This script defines square segments of Greenland's coast that contain
% any densified terminus position data. 
% 
% Chad Greene, NASA/JPL, August 2022. 

%% Load densified terminus picks
% Load data that was created by terminus_data_densifier.m. 

T = load('/Users/cgreene/Documents/GitHub/greenland-coastlines/data/terminus_data_densified.mat'); 

xi = double(T.x); 
yi = double(T.y); 
ti = double(T.t); 
[yeari,~,~] = datevec(ti); 

%% Create a grid of where teminus data exist
% Use gridbin (on GitHub) to create a Greenland-wide mask where true cells
% contain terminus data, and false everywhere else. This step is much more
% efficient than methods like kmeans clustering. Larger grid resolution is
% more efficient (in this script), but creates larger polygons that may
% cover too large of a geographic area. In my tests, 5 km looks about right.
% 
% After creating the binary grid G, use bwlabel to find the connected
% segments of true grid cells in G. This ensures any terminus position that
% straddles multiple grid cells will be counted in the same segment. 

res = 5e3; % (m) Grid resolution 

% Define the grid coordinates: 
x = (min(xi)-res/2):res:(max(xi)+res/2);
y = (max(yi)+res/2):-res:(min(yi)-res/2);

% A grid of where we have data:
G = gridbin(xi,yi,true(size(xi)),x,y,@any); 

% Label the connected grid cells: 
[L,n] = bwlabel(G,4); 

figure
h=imagescn(x,y,L);
h.AlphaData = G; 
axis image off 

%% Define boundaries of each cluster of terminus position data 
% Using the spatial extents of the data within each cluster of terminus
% positions, define a rectangle for each cluster that spans the data plus a
% specified buffer. 

buf = 5e3; % (m) buffer on all sides of measured position data 

% Determine which pixel cluster each terminus position point belongs to:
Li = interp2(x,y,L,xi,yi,'nearest'); 

% Preallocate:
xmin = nan(n,1); 
xmax = xmin; 
ymin = xmin; 
ymax = xmin; 

% Define the extents of each pixel cluster:
for k = 1:n
%     % The total extents of all data in this cluster: 
    xmin(k) = min(xi(Li==k)) - buf; 
    xmax(k) = max(xi(Li==k)) + buf; 
    ymin(k) = min(yi(Li==k)) - buf; 
    ymax(k) = max(yi(Li==k)) + buf; 
end

% Center coordinates: 
xc = (xmin+xmax)/2; 
yc = (ymin+ymax)/2; 

% Pad to square aspect ratio:
dx = xmax-xmin; 
dy = ymax-ymin; 

widen = dy>dx; % indices of rectangles to make wider 
xmin(widen) = xc(widen) - dy(widen)/2; 
xmax(widen) = xc(widen) + dy(widen)/2; 

heighten = dx>dy; % indices of rectangles to make taller
ymin(heighten) = yc(heighten) - dx(heighten)/2; 
ymax(heighten) = yc(heighten) + dx(heighten)/2; 

% Convert to polyshape: 
for k = 1:n
    P(k) = polyshape([xmin(k) xmax(k) xmax(k) xmin(k)],[ymin(k) ymin(k) ymax(k) ymax(k)]); 
end


%%

figure('pos',[50 50 560 840])
plot(P)
axis image off
hold on
fastscatter(xi,yi,yeari,'markersize',1)
grid off
cb = colorbar; 
caxis([1985 2020])
modismog('contrast','white')
cmocean thermal

%export_fig('/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/densified_data_segmenter.png','-r600','-p0.01');
axis([-256493     -159058    -2320514    -2156606])
%export_fig('/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/densified_data_segmenter_zoom.png','-r600','-p0.01');

%%

for k = 1:n
    S(k).Geometry = 'Polygon'; 
    S(k).X = [xmin(k) xmax(k) xmax(k) xmin(k) xmin(k)];
    S(k).Y = [ymin(k) ymin(k) ymax(k) ymax(k) ymin(k)]; 
    S(k).projcrs = 3413;
end

%shapewrite(S,'/Users/cgreene/Documents/GitHub/greenland-coastlines/data/calving_front_polygons.shp')
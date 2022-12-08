

%% Load data 

T = load('terminus_data_densified_2022-11-14.mat'); 

T.x = double(T.x); 
T.y = double(T.y); 
T.t = double(T.t); 

[T.t,ind] = sort(T.t); 
T.x = T.x(ind); 
T.y = T.y(ind); 
T.p = T.p(ind); 

%%

axl=[-386943.36    -314366.74   -1651512.22   -1500578.16; % tall
   -336984.22    -293337.91   -1862499.52   -1754878.89;
   -280404.35    -223945.45   -1960522.23   -1900052.21; % single
   -250964.42    -165793.80   -2289442.38   -2164237.40; %jako
   -233137.03    -157355.26   -3156088.86   -3052689.24; % nothing sw
    71547.79     143198.86   -3236549.87   -3128168.11; % not much se
   80751.82     139030.22   -3058624.53   -2983725.28; % nice
   80751.82     139030.22   -3058624.53   -2983725.28; % also nice
   176543.34     282329.95   -2759348.15   -2642329.44; % ok if needed
   283378.91     440430.42   -2609628.46   -2535290.73; %wide
   444853.57     735013.89   -2337715.08   -2185241.26; % could be ok
   -445406.68    -232449.19   -1107868.93    -894911.44; % petermann
   -630413.91    -439832.25   -1376756.58   -1225848.86;
   -584408.33    -426314.60   -1432997.75   -1380299.84];

%%

cax = [datenum(1985,1,1) datenum(2022,1,1)]; 

figure('color','k')

% top left
ax(1) = axes('position',[0 1/3 1/4 2/3]);
axis(axl(1,:))
axis([mean(xlim)+diff(ylim)*[-1 1]/4 mean(ylim)+diff(ylim)*[-1 1]/2])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([mean(xlim)+diff(ylim)*[-1 1]/4 mean(ylim)+diff(ylim)*[-1 1]/2])
txtcb=textcolorbar(1985:2022,'location','ne','colormap',cmocean('thermal'),'fontsize',5);
[hsb,hsbt] = scalebarpsn('len',10,'location','sw','color','w','fontsize',5); 
hsb.LineWidth = 1; 
hsb.YData = hsb.YData-3e3; 
hsbt.Position(2) = hsbt.Position(2)-3e3; 

ntitle('a','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)
txtcb.BackgroundColor=[.5 .5 .5 .5];
%txtcb.BackgroundColor=[1 1 1 .5];
txtcb.Margin=.0001;

% m(1) = mapzoompsn_coast('ne','frame','off'); 
% m(1).Position(1) = 0.155;

% bottom left
ax(2) = axes('position',[0 0 3/4 1/3]);
axis(axl(14,:)+5e3*[-1 1 -1 1])
%axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/4])
axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/6])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/6])
[hsb,hsbt] = scalebarpsn('len',10,'location','sw','color','w','fontsize',5); 
hsb.LineWidth = 1; 
ntitle('g','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)
%m(2) = mapzoompsn_coast('ne','frame','off','insetsize',0.5); 


% Petermann
ax(3) = axes('position',[1/4 2/3 1/4 1/3]);
axis(axl(12,:))
axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/2])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/2])
[hsb,hsbt] = scalebarpsn('len',20,'location','se','color','k','fontsize',5); 
hsb.LineWidth = 1; 
ntitle('b','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)
%m(3) = mapzoompsn_coast('ne','frame','off','insetsize',0.5); 
%
% fuse
ax(4) = axes('position',[2/4 2/3 1/4 1/3]);
% axis(axl(3,:))
% axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/2])
axis([ -280065.80    -226306.42   -1957682.64   -1903923.26])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([ -280065.80    -226306.42   -1957682.64   -1903923.26])
[hsb,hsbt] = scalebarpsn('len',10,'location','se','color','k','fontsize',5); 
hsb.LineWidth = 1; 
ntitle('c','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)


% bottom right 
ax(5) = axes('position',[3/4 0 1/4 1/3]);
% axis(axl(9,:))
% axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/2])
%axis([75959.19     122130.97   -3045816.13   -2999644.35])
%axis([ 172867.03     275317.78   -2757798.14   -2655299.10])
axis([ 172867.03     275317.78 mean([-2757798.14   -2655299.10])+[-1 1]*diff([172867.03     275317.78])/2])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
%axis([75959.19     122130.97   -3045816.13   -2999644.35])
%axis([ 172867.03     275317.78   -2757798.14   -2655299.10])
axis([ 172867.03     275317.78 mean([-2757798.14   -2655299.10])+[-1 1]*diff([172867.03     275317.78])/2])
[hsb,hsbt] = scalebarpsn('len',10,'location','se','color','w','fontsize',5); 
hsb.LineWidth = 1; 
ntitle('h','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)

% top right
ax(6) = axes('position',[3/4 2/3 1/4 1/3]);
%axis(axl(7,:))
% axis([76208.30     124955.76   -3047533.83   -2998786.38])
% axis([mean(xlim)+diff(xlim)*[-1 1]/2 mean(ylim)+diff(xlim)*[-1 1]/2])
axis([75959.19     122130.97   -3045816.13   -2999644.35])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([75959.19     122130.97   -3045816.13   -2999644.35])
[hsb,hsbt] = scalebarpsn('len',10,'location','sw','color','k','fontsize',5); 
hsb.LineWidth = 1; 
ntitle('d','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)


% middle right
% %ax(7) = axes('position',[3/4 1/3 1/4 1/3]);
% ax(7) = axes('position',[2/4+1/8 1/3 3/8 1/3]);
% %axis(axl(10,:))
% axis([291355.02     419003.86   -2627609.26   -2499960.42])
% axis([mean(xlim)+diff(xlim)*[-1 1] mean(ylim)+diff(xlim)*[-1 1]/2])
% drawnow
% plotstuff(xlim,ylim,xc,yc,col)
% axis([mean(xlim)+diff(xlim)*[-1 1] mean(ylim)+diff(xlim)*[-1 1]/2])

%ax(7) = axes('position',[3/4 1/3 1/4 1/3]);
ax(7) = axes('position',[2/4+1/8 1/3 3/8 1/3]);
xi=359065;
yi=-2563317; 
dx = 67e3;
xi=356512;
yi=-2563955; 
dy = dx/1.5; 
axis([xi-dx xi+dx yi-dy yi+dy])

drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([xi-dx xi+dx yi-dy yi+dy])
%axis([291283.95     412416.92   -2609779.22   -2529023.91])
[hsb,hsbt] = scalebarpsn('len',10,'location','sw','color','k','fontsize',5); 
hsb.LineWidth = 1; 
ntitle('f','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)


ax(8) = axes('position',[1/4 1/3 3/8 1/3],'color','none');
xi=359065;
dx = 38e3/2;
xi=-192597-1e3;
yi=-2270908-1e3; 
dy = dx/1.5; 
axis([xi-dx xi+dx yi-dy yi+dy])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([xi-dx xi+dx yi-dy yi+dy])
%axis([291283.95     412416.92   -2609779.22   -2529023.91])
[hsb,hsbt] = scalebarpsn('len',5,'location','se','color','k','fontsize',5); 
hsb.LineWidth = 1; 

ntitle('e','location','northwest','fontweight','bold','fontsize',7,'color','w','backgroundcolor','k','margin',.01)

% export_fig('/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/greenland_coastline_summary_2022-11-17.jpg','-r600','-opengl','-nocrop','-p0.002')

%%
if false 

figure('pos',[47   375   729   572])

xi=359065;
dx = 38e3/2;
xi=-192597-1e3;
yi=-2270908-1e3; 
dy = dx/1.5; 
axis([xi-dx xi+dx yi-dy yi+dy])
drawnow
plotstuff(xlim,ylim,T.x,T.y,T.t,cax)
axis([xi-dx xi+dx yi-dy yi+dy])
%axis([291283.95     412416.92   -2609779.22   -2529023.91])
[hsb,hsbt] = scalebarpsn('len',5,'location','se','color','k','fontsize',8); 
hsb.LineWidth = 1; 
set(gcf,'color','k')
txtcb=textcolorbar(1985:2021,'location','e','colormap',cmocean('thermal'),'fontsize',8);

end

%% SUBFUNCTIONS

function plotstuff(xl,yl,xc,yc,tc,cax,p) 


   [Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xl,yl,1000);
   image(xx,yy,Itmp(:,:,1:3));
   daspect([1 1 1])
   hold on

   if nargin>6
       for k = 1:max(p)
           ind = xc>=xl(1) & xc<=xl(2) & yc>=yl(1) & yc<=yl(2) & p==k; 
           fastscatter(xc(ind),yc(ind),tc(ind),'markersize',0.1)
       end
   else
       ind = xc>=xl(1) & xc<=xl(2) & yc>=yl(1) & yc<=yl(2);
       fastscatter(xc(ind),yc(ind),tc(ind),'markersize',0.1)
   end

   cmocean thermal

   caxis(cax)
   axis xy 
   label_greenland_glaciers('fontsize',5,'color',0.1*[1 1 1],'shadow','abbreviate')
   box on
   set(gca,'xtick',[],'ytick',[],'linewidth',2,'xcolor','k','ycolor','k')

end





















function h = mapzoompsn_coast(varargin)
% mapzoompsn zooms a north polar stereographic map to a specified location and extent
% and/or places an inset map of Greenland for spatial context. 
% 
%% Citing Antarctic Mapping Tools
% This function was adapted from Antarctic Mapping Tools for Matlab (AMT). If it's useful for you,
% please cite our paper: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
% @article{amt,
%   title={{Antarctic Mapping Tools for \textsc{Matlab}}},
%   author={Greene, Chad A and Gwyther, David E and Blankenship, Donald D},
%   journal={Computers \& Geosciences},
%   year={2017},
%   volume={104},
%   pages={151--157},
%   publisher={Elsevier}, 
%   doi={10.1016/j.cageo.2016.08.003}, 
%   url={http://www.sciencedirect.com/science/article/pii/S0098300416302163}
% }
%   
%% Syntax
% 
% mapzoompsn
% mapzoompsn(lat,lon) 
% mapzoompsn(x,y)
% mapzoompsn(...,'size',mapsizekm)
% mapzoompsn(...,InsetLocation)
% mapzoompsn(...,'insetsize',sizefraction)
% mapzoompsn(...,'frame','off')
% mapzoompsn(...,'km')
% h = mapzoompsn(...) 
% 
%% Description 
% 
% mapzoompsn(lat,lon) centers a 500 km wide map about the georeferenced 
% location given by lat, lon. 
% 
% mapzoompsn(x,y) centers a 500 km wide map about the polar stereographic 
% eastings and northings x and y. 
% 
% mapzoompsn(...,'size',mapsizekm) specifies size of the map in kilometers given 
% mapsizekm, which can be a scalar to create a square map or a two-element array
% for a rectangular map in the form [mapwidthkm mapheightkm], where mapwidthkm and
% mapheightkm are the dimensions of the map in kilometers. 
%
% mapzoompsn(...,InsetLocation) creates an inset map at the location InsetLocation,
% which can be 
%           'southeast' or 'se'  lower right corner 
%           'northwest' or 'nw'  upper left corner
%           'northeast' or 'ne'  upper right corner
%           'southwest' or 'sw'  lower left corner
%
% mapzoompsn(...,'insetsize',sizefraction) specifies size of the inset as a
% fraction of the width of the current map. Default sizefraction is 0.25. 
%
% mapzoompsn(...,'frame','off') removes frame from the inset. 
%
% mapzoompsn(...,'km') is for plots in polar stereographic kilometers rather than the default meters.
% 
% h = mapzoompsn(...) returns a handle h of inset map axes. 
% 
%% Example 1 
% Zoom in on Petermann Glacier like this: 
% 
%   greenland
%   mapzoompsn(80.75,-65.75,'ne')
%
%% Author Info 
% This function and supporting documentation were written by Chad A. Greene of the 
% University of Texas at Austin's Institute for Geophysics (UTIG), June 2017. 
% Feel free to contact me if you have any questions or comments. 
% http://www.chadagreene.com
% 
% See also scarloc, scarlabel, scarclick, and scalebarpsn. 

%% Set defaults: 

inset = false; 
insetsize = 0.25; 
frameon = true; 
location = 'northeast'; 
usekm = false; 
if nargin==0 
   UseCurrentExtents = true; 
else
   UseCurrentExtents = false; 
   mapsize = [500 500]; % sets default map size to 500 km by 500 km
end

%% Parse inputs: 

% Inset location: 
tmp = strcmpi(varargin,'southwest')|strcmpi(varargin,'northwest')|...
      strcmpi(varargin,'southeast')|strcmpi(varargin,'northeast')|...
      strcmpi(varargin,'sw')|strcmpi(varargin,'nw')|...
      strcmpi(varargin,'se')|strcmpi(varargin,'ne'); 
if any(tmp)
   inset = true; 
   location = varargin{tmp}; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
end

% Check for inset size declaration: 
tmp = strcmpi(varargin,'insetsize'); 
if any(tmp) 
   inset = true; 
   insetsize = varargin{find(tmp)+1}; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
end
   
% Check for frame declaration: 
tmp = strcmpi(varargin,'frame');
if any(tmp)
   inset = true; 
   if strcmpi(varargin{find(tmp)+1},'off')||strcmpi(varargin{find(tmp)+1},'none');
      frameon = false; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
   end
end

% Map width: 
tmp = strcmpi(varargin,'size')|strcmpi(varargin,'mapsize')|strcmpi(varargin,'mapwidth')|strcmpi(varargin,'mapwidthkm')|strcmpi(varargin,'width'); 
if any(tmp)
   mapsize = varargin{find(tmp)+1}; 
   assert(isnumeric(mapsize)==1,'Map size must be numeric.'); 
   if isscalar(mapsize) 
      mapsize = [mapsize mapsize]; 
   end
   assert(numel(mapsize)==2,'Map size must be a one- or two-element numeric value.') 
end

% Polar stereographic kilometers or meters? 
tmp = strcmpi(varargin,'km'); 
if any(tmp) 
   usekm = true; 
   if tmp(1)
      UseCurrentExtents = true; 
   end
end
   

% Center location declaration: 
if ~UseCurrentExtents

  % User has entered location by coordinates: 
  if islatlon(varargin{1},varargin{2})
     [xc,yc] = ll2psn(varargin{1},varargin{2}); 
  else
     xc = varargin{1}; 
     yc = varargin{2}; 
  end
   
end

%% Set axes of map: 

gcah = gca; % handle of initial plot

if UseCurrentExtents
   ax = axis; 
else
   axis equal xy
   ax = [xc-mapsize(1)*500 xc+mapsize(1)*500 yc-mapsize(2)*500 yc+mapsize(2)*500]; 
   axis(ax); 
end

% Define x,y coordinates of axis perimeter: 
axx = [ax(1) ax(2) ax(2) ax(1) ax(1)];
axy = [ax(3) ax(3) ax(4) ax(4) ax(3)]; 
   
%% Place an inset map: 

if inset
      
    gp = plotboxpos(gca); 
    
    insetwidth = insetsize*gp(3); 
    insetheight = insetsize*gp(4); % just changed this to 4 (was 3 for a year or so? ) 
    
    switch lower(location)
        case {'southwest','sw'}
            insetx = gp(1);
            insety = gp(2);   
            
        case {'northeast','ne'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'northwest','nw'}
            insetx = gp(1); 
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'southeast','se'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2);            
            
        otherwise 
            error('Unrecognized inset location.')
    end
    

   % Create new set of axes for inset map: 
   h = axes('position',[insetx insety insetwidth insetheight],'tag','insetmap');
   hold on
   
  
   % Plot greenland: 
   greenland('patch','facecolor',0.75*[1 1 1],'linewidth',0.2)
   
   % Plot red box:
   if usekm
      plot(axx*1000,axy*1000,'r-','linewidth',1); 
   else
      plot(axx,axy,'r-','linewidth',1); 
   end
      
   axis equal tight
   
   
   % Set final dimensions after plotting the inset: 
   gpinset = plotboxpos(gca); 
   insetwidth = gpinset(3); 
   insetheight = gpinset(4); % just changed this to 4 (was 3 for a year or so? ) 
    
    switch lower(location)
        case {'southwest','sw'}
            insetx = gp(1);
            insety = gp(2);   
            
        case {'northeast','ne'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'northwest','nw'}
            insetx = gp(1); 
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'southeast','se'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2);            
            
        otherwise 
            error('Unrecognized inset location.')
    end
   
   
   % Format inset axes: 
   set(gca,'xtick',[],'ytick',[],'position',[insetx insety insetwidth insetheight])
   if frameon
      box on
   else
      axis off
   end
   
  
   % Make the original map axes current: 
   axes(gcah); 
   
   % Ensure inset map is on top of the stack: 
   uistack(gcah,'down');
   
   
  % Clean up: 
  if nargout==0 
     clear h
  end
end


end

%% Kelly Kearney's plotboxpos function: 

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);
end        
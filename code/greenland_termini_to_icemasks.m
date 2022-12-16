% This script turns terminus_data_densified data into a time-evolving
% monthly ice mask for Greenland. This script has been adapted from
% termpicks2mask_v2.m.

devon = false; % Set to true if running the whole ice sheet on Devon. 

% datenumbers: The first data in 1972 is Sept 6
t = datenum(datetime(1972,9,15):calmonths(1):datetime(2022,2,15));

readme = ['Written by ',mfilename,'.m.']

if devon    
   addpath(genpath('/mnt/devon2-r1/devon0/cgreene/greenland_masking'))
   cd('/usr/local/matlab-9.9/toolbox/matlab/specgraph/private')
   fn_term = '/mnt/devon2-r1/devon0/cgreene/greenland_masking/terminus_data_densified_2022-11-14.mat';
   fn_initial = ['/mnt/devon2-r1/devon0/cgreene/greenland_masking/greenland_icemask_initial_',datestr(now,'yyyy-mm-dd'),'.mat']; 
   fn_backfill = ['/mnt/devon2-r1/devon0/cgreene/greenland_masking/greenland_icemask_backfill_',datestr(now,'yyyy-mm-dd'),'.mat'];
   fn_prefinal = ['/mnt/devon2-r1/devon0/cgreene/greenland_masking/greenland_icemask_prefinal_',datestr(now,'yyyy-mm-dd'),'.mat'];
   fn_final = ['/mnt/devon2-r1/devon0/cgreene/greenland_masking/greenland_icemask_final_',datestr(now,'yyyy-mm-dd'),'.mat'];
   fn_clean = ['/mnt/devon2-r1/devon0/cgreene/greenland_masking/greenland_icemask_clean_',datestr(now,'yyyy-mm-dd'),'.mat'];
   newfilename = ['/mnt/devon2-r1/devon0/cgreene/greenland_masking/greenland_monthly_ice_masks_',datestr(now,'yyyy-mm-dd'),'.nc']; 

else
   cd('/Applications/MATLAB_R2022a.app/toolbox/matlab/specgraph/private')
   fn_term = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/terminus_data_densified_2022-11-14.mat';
   fn_initial = ['/Users/cgreene/Documents/MATLAB/greenland_icemask_initial_',datestr(now,'yyyy-mm-dd'),'.mat']; 
   fn_backfill = ['/Users/cgreene/Documents/MATLAB/greenland_icemask_backfill_',datestr(now,'yyyy-mm-dd'),'.mat'];
   fn_prefinal = ['/Users/cgreene/Documents/MATLAB/greenland_icemask_prefinal_',datestr(now,'yyyy-mm-dd'),'.mat'];
   fn_final = ['/Users/cgreene/Documents/MATLAB/greenland_icemask_final_',datestr(now,'yyyy-mm-dd'),'.mat'];
   fn_clean = ['/Users/cgreene/Documents/MATLAB/greenland_icemask_clean_',datestr(now,'yyyy-mm-dd'),'.mat'];
   newfilename = ['/Users/cgreene/Documents/MATLAB/greenland_monthly_ice_masks_',datestr(now,'yyyy-mm-dd'),'.nc']; 
end

%% Load terminus positions: 
% These terminus positions have already been compiled from multiple
% datasets and densified by terminius_data_densifier.m: 

T = load(fn_term); 

T.x = double(T.x); 
T.y = double(T.y); 
T.t = double(T.t); 

% Delete any potential duplicates: 
[~,ind] = unique([T.t' T.x' T.y' double(T.p)'],'rows'); 
T.t = T.t(ind); 
T.x = T.x(ind); 
T.y = T.y(ind); 
T.p = T.p(ind); 

% Delete anything long before the first data cube slice: 
ind = T.t>=(t(1)-31); 
T.t = T.t(ind); 
T.x = T.x(ind); 
T.y = T.y(ind); 
T.p = T.p(ind); 

clear ind 

%% Load gridded data 

filename = 'greenland_extruded_velocity_and_thickness_2022-12-15.nc'; 

x = double(ncread(filename,'x')); 
y = double(ncread(filename,'y')); 

vx = double(permute(ncread(filename,'vx'),[2 1])); 
vy = double(permute(ncread(filename,'vy'),[2 1])); 
rock = permute(ncread(filename,'v_source'),[2 1])==0; 
th = double(permute(ncread(filename,'thickness'),[2 1])); 

[X,Y] = meshgrid(x,y); 
bed = bedmachine_interp('bed',X,Y,'greenland'); 
marine = bed<0;

fn = 'GimpOceanMask_90m_2015_v1.2.tif'; 
[mx,my]=pixcenters(geotiffinfo(fn));
tmp = filt2(double(imread(fn)),90,250,'lp'); % lowpass filter before interpolation to anti-alias 
ocean2015 = interp2(mx,my,tmp,X,Y,'linear',1)>=0.5;
ocean2015(rock) = false; 
ocean2015(bed>0) = false; 

land2015 = imfill(~ocean2015,8,'holes'); % any accidental holes created when interpolating 
clear fn mx my ocean2015 tmp filename bed

%% Trim the grid for testing
% This is only to reduce grid size for quick tests. 

if devon
    trimtest = false; 
else
    trimtest = 'jako'; % can be 'jako', 'nw', 'slim', or false 
end

switch trimtest
   
   case 'jako' % Jakobshavn, same extents as extruded velocity masking figure.
      cols = x>=-244567.60  & x<=-173694.31; 
      rows = y>=-2334333.46 & y<=-2210298.49;
      
   case 'nw' % The northwest, downsampled significantly 
      rows = 3000:4:8000; 
      cols = 500:4:4000; 
      
   case 'slim' % The entire landmass of greenland and not much more. This option trims about 11% of the total number of pixels.  
      cols = x>=-637310  & x<=849704; 
      rows = y>=-3355603 & y<=-664846;
      
   case false 
      % do nothing
      
   otherwise
      error('unrecognized trimtest region')
end

if trimtest % will be true for 'jako' or 'nw' 

   x = x(cols); 
   y = y(rows); 
   vx = vx(rows,cols); 
   vy = vy(rows,cols); 
   th = th(rows,cols); 
   rock = rock(rows,cols); 
   marine = marine(rows,cols); 
   land2015 = land2015(rows,cols); 
   X = X(rows,cols); 
   Y = Y(rows,cols); 
   
   keep = T.x>min(x) & T.x<max(x) & T.y>min(y) & T.y<max(y); 
   T.x = T.x(keep); 
   T.y = T.y(keep); 
   T.t = T.t(keep); 
   T.p = T.p(keep); 
   
   clear rows cols keep 
end

yf = flipud(y); 
vxf = flipud(vx); 
vyf = flipud(vy); 

%% Manual edits 

% Starting in Jun 2003 (k=1052) an embayment near Jakobshavn clears out, but not all terminus positions capture that, so we intervene manually.   
jakotrim = inpolygon(X,Y,...
   [-202692     -202321     -192994     -192013     -191324     -191430],...
   [-2273212    -2280208    -2280208    -2277319    -2275756    -2274060]); 

%% Prepare land cube

% Preallocate land cube for this year with the correct number of time slices:: 
land = false(size(X,1),size(X,2),numel(t));

%% Create example image 
% This code section creates a cartoon of mask adjustments. 

SaveAdjustmentAnimations = false; % Don't do this unless you want to overwrite current files.  

if SaveAdjustmentAnimations

    step = 0.2; % Number of streamline steps per grid cell.  
    reach = 50e3; % meters of reach in each direction from measured termini
    Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 
       
    % March forward in time starting aug 2015, then march backward 
    k0 = find(t==datenum(datetime('aug 15, 2015')));
    
    k=k0; 
    
    % This month's terminus datapoints: 
    indp = T.t>=(t(k)-30) & T.t<=t(k); % Terminus points within the past 30 days
    indf = T.t<=(t(k)+30) & T.t>=t(k); % Terminus points within the next 30 days
    
    % Get terminus positions taken within 30 days in the past:  
    cxp = T.x(indp);
    cyp = T.y(indp);
    cpp = T.p(indp); 
    dtp = T.t(indp) - t(k); 
    
    % Account for <=30 days of advection:
    dxp = interp2(x,y,vx,cxp,cyp) .* dtp/365.25;
    dyp = interp2(x,y,vy,cxp,cyp) .* dtp/365.25;
    cxp = cxp - dxp; 
    cyp = cyp - dyp; 
    
    % Get terminus positions taken less than 30 days in the future:  
    cxf = T.x(indf);
    cyf = T.y(indf);
    cpf = T.p(indf); 
    dtf = T.t(indf) - t(k); 
    
    % Account for <=30 days of advection:
    dxf = interp2(x,y,vx,cxf,cyf) .* dtf/365.25;
    dyf = interp2(x,y,vy,cxf,cyf) .* dtf/365.25;
    cxf = cxf - dxf; 
    cyf = cyf - dyf; 
        
    kp=4; 
      
    % Indices of advected past points that fall within unconstrained grid cells:
    indkpf = cpf==kp & cxf>x(1) & cxf<x(end) & cyf<y(1) & cyf>y(end); 
    
    % Upstream from future positions advected to today: 
    st = stream2(x,y,-vx,-vy,cxf(indkpf),cyf(indkpf),[step Npts]);
    us = cell2mat(st(:)); % stream paths upstream of the terminus
    if ~isempty(us)
        isf = isfinite(us(:,1)); 
        Us = gridbin(us(isf,1),us(isf,2),true(size(us(isf,1))),x,y,@any); % gridded cells upstream of the terminus 
        Us = imfill(bwmorph(Us,'close'),4,'holes'); 
    else
        Us = false(size(X)); 
    end
    
    % Indices of advected past points that fall within unconstrained grid cells:
    indkpp = cpp==kp & interp2(x,y,~rock,cxp,cyp,'nearest') & cxp>x(1) & cxp<x(end) & cyp<y(1) & cyp>y(end); 
    
    % Downstream from past positions advected to today: 
    st = stream2(x,y,vx,vy,cxp(indkpp),cyp(indkpp),[step Npts]);
    ds = cell2mat(st(:)); % stream paths from terminus 
    if ~isempty(ds)
        isf = isfinite(ds(:,1)); 
        Ds = gridbin(ds(isf,1),ds(isf,2),true(size(ds(isf,1))),x,y,@any); % Grid cells that are downstream of the terminus 
        Ds = imfill(bwmorph(Ds,'close'),4,'holes'); 
    else
        Ds = false(size(X)); 
    end
      
    
    FillRegion = Us & ~rock; 
    CarveRegion = Ds & ~rock; 
    
    tmp = land2015; 
    tmp(CarveRegion) = false; 
    tmp(FillRegion) = true; 
    tmp(rock) = true; 
    tmp = imfill(tmp,4,'holes'); 
    CarveRegion(tmp) = false; 
    
    %
    figure('pos',[20 20 270 205])
    h1=imagescn(x,y,land2015 | rock); % unadjusted  
    h1.AlphaData = ~rock; 
    axis([-197879     max(x) -2281327    -2261633])
    daspect([1 1 1]) 
    axis  off
    hold on
    cmocean ice 
    
    [Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xlim,ylim,1000);
    hold on
    hh=image(xx,yy,Itmp(:,:,1:3));
    uistack(hh,'bottom'); 
    set(gca,'pos',[0 0 1 1])
    
    hpds = plot(cxp(indkpp),cyp(indkpp),'.','color',hex2rgb('#F6F930'),'markersize',2); 
    hpus = plot(cxf(indkpf),cyf(indkpf),'.','color',hex2rgb('#0496ff'),'markersize',2); 
    hds = maskoverlay(x,y,CarveRegion,'color',hpds.Color,'alpha',0.5);
    hus = maskoverlay(x,y,FillRegion,'color',hpus.Color,'alpha',0.5);
    uistack(hpds,'top'); 
    uistack(hpus,'top');
    
    txt1 = text(-189796,-2272091,'Carve region','horiz','center','vert','middle',...
        'fontsize',8,'color',hpds.Color,'fontweight','bold'); 
    
    txt2 = text(-180599,-2270202.23,'Fill region','horiz','center','vert','middle',...
        'fontsize',8,'color',hpus.Color,'fontweight','bold'); 
    
    txta = text(-182158,-2272183,'a','horiz','left','vert','bot','fontsize',6);% fill 
    txtb = text(-182788,-2277963,'b','horiz','left','vert','bot','fontsize',6); % carve 
    txtc = text(-180673,-2277534,'c','horiz','center','vert','mid','fontsize',6); % no change 
    txtd = text(-180165,-2274773,'d','horiz','left','vert','bot','fontsize',6); % overlap
    
    [sb1,sb2] = scalebarpsn('location','se','fontsize',6,'length',1);
    sb1.LineWidth = 1; 
    
    export_fig('/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/carve_fill_cartoon.jpg','-r600')
    
    txt_title = ntitle('Prior','fontweight','bold','fontsize',8); 
    gif('/Users/cgreene/Documents/GitHub/greenland-coastlines/animations/carve_fill_cartoon.gif','resolution',900,'frame',gca,'delaytime',1)
    
    txt_title.String = 'Adjusted'; 
    h1.CData = tmp; 
    gif
    
    clear sb* txt* h1 hp*
end

%% Create initial cube

step = 0.2; % Number of streamline steps per grid cell.  
reach = 50e3; % meters of reach in each direction from measured termini
Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 

disp([datestr(now),': Starting initial land cube.'])
   
% March forward in time starting aug 2015, then march backward 
k0 = find(t==datenum(datetime('aug 15, 2015')));

sometimes_ice = false(size(rock)); 

for k = [k0:length(t) (k0-1):-1:1]

   % Determine which mask to adjust: 
   % Start with the most recently-adjusted mask, which we will adjust again to make today's mask: 
   switch sign(k-k0)
      case 0 
         tmp = land2015; 
         clear land2015
      case 1
         % We're moving forward in time, so start with last month's mask 
         tmp = land(:,:,k-1); 
      case -1 
         % We're moving backward in time, so start with next month's mask
         tmp = land(:,:,k+1); 
   end
    
    % This month's terminus datapoints: 
    ind = abs(T.t-t(k))<=16;

    if any(ind)
       
       % Get terminus positions taken within 30 days in the past:  
       cx = T.x(ind);
       cy = T.y(ind);
       cp = T.p(ind); 
       dt = T.t(ind) - t(k); 
       
       % Account for <=30 days of advection:

       % Building griddedInterpolant once is much faster than 
       if k==k0
          Fx = griddedInterpolant({yf,x},vxf,'linear','none'); 
          Fy = griddedInterpolant({yf,x},vyf,'linear','none'); 
       end

       dx = Fx(cy,cx) .* dt/365.25;
       dy = Fy(cy,cx) .* dt/365.25;
       
       cx = cx - dx; 
       cy = cy - dy; 
       
       % Loop through all priorities of data, starting with the lousiest: 
       for kp = 1:max(T.p)
   
           % Indices of advected past points that fall within unconstrained grid cells:
           %indkp = cp==kp &  interp2(x,y,~rock,cx,cy,'nearest') & cx>x(1) & cx<x(end) & cy<y(1) & cy>y(end); 
           indkp = cp==kp & cx>x(1) & cx<x(end) & cy<y(1) & cy>y(end); 
   

           % Upstream from future positions advected to today: 
           if any(indkp)
              [xus,yus] = stream2_fast(x,yf,-vxf,-vyf,cx(indkp),cy(indkp),[step Npts]);
           else 
              xus = []; 
           end

           if ~isempty(xus)
               isf = isfinite(xus); 
               Us = gridbin(xus(isf),yus(isf),true(size(xus(isf))),x,y,@any); % gridded cells upstream of the terminus 
               Us = imfill(bwmorph(Us,'close'),4,'holes'); 
           else
               Us = false(size(X)); 
           end 


           % Downstream from past positions advected to today: 
           if any(indkp)
              [xds,yds] = stream2_fast(x,yf,vxf,vyf,cx(indkp),cy(indkp),[step Npts]);
           else 
              xds = []; 
           end
           if ~isempty(xds)
               isf = isfinite(xds); 
               Ds = gridbin(xds(isf),yds(isf),true(size(xds(isf))),x,y,@any); % Grid cells that are downstream of the terminus 
               Ds = imfill(bwmorph(Ds,'close'),4,'holes'); 
           else
               Ds = false(size(X)); 
           end
             
           FillRegion = Us & ~rock; 
           CarveRegion = Ds & ~rock; 
           sometimes_ice(FillRegion) = true; 
           
           tmp(FillRegion) = true; 
           tmp(CarveRegion) = false; % overwrites any potential FillRegion
           tmp(imdilate(tmp,strel('disk',3)) & Us & Ds & ~tmp) = true; % Defines the ice edge as true 
           tmp(rock) = true; 
           tmp = imfill(tmp,8,'holes'); 
           
       end
    end

    land(:,:,k) = tmp; 
    disp([datestr(now),': Initial cube finished month ',num2str(k),' of ',num2str(length(t)),' with k0=',num2str(k0)])

end


always_land = all(land,3); 
never_land = ~any(land,3) & ~sometimes_ice;

always_land(rock) = true; 
never_land(rock) = false; 

sometimes_ice = ~always_land & ~never_land; 

% Fix jagged edges: 
for k = 1:size(land,3)
   tmp = land(:,:,k); 
   tmp = bwmorph(bwmorph(tmp,'close'),'open'); 
   tmp(always_land) = true; 
   tmp(never_land) = false;
   land(:,:,k) = tmp; 
end

land = remove_icebergs(land,marine,always_land,never_land);

disp([datestr(now),': Done masking the initial cube.'])

landsum = squeeze(sum(sum(land))); 

clear CarveRegion FillRegion reach cx cy cp dx dy dt cpf cpp st cxf cxp cyf cyp us ds Ds dtf dtp dxf dxp dyf dyp Us ind* isf Itmp k k0 kp tmp trimtest st  

% Write the data in case of crash: 
save(fn_initial,'-v7.3')

%% Backfill 1 month
% Work backwards. Calculate expected displacement from one month to the
% next. If, after accounting for advection, next months's mask says this
% month's grid cells must have been ice, fill this month's grid cells with
% ice. We're only doing one month in this particular cell, because here
% we're accounting for variations in dt (10%) due to different lengths of
% months. 

vtol = 1.0; % velocity tolerance scale factor. 1.05 allows 5% wiggle room for velociy. 1.0 strictly enforces velocity 

land_backfill = land; 
dt = diff(t); 

% Daily displacement: 
dX = vtol*interp2(x,y,vx,X,Y)/365.25; 
dY = vtol*interp2(x,y,vy,X,Y)/365.25; 

for k = (numel(t)-1):-1:1
   
   % Find out if it will be ice next month by interpolating next month's mask at today's grid point positions plus dt days of displacement: 
   %will_be_ice_old = interp2(x,y,single(land_backfill(:,:,k+1)),X+dX*dt(k),Y+dY*dt(k))==1;

   % Building griddedInterpolant once is much faster than 
   if k==(numel(t)-1)
      F = griddedInterpolant({flipud(y),x},single(flipud(land_backfill(:,:,k+1))),'linear','none'); 
   else
      F.Values = single(flipud(land_backfill(:,:,k+1)));
   end
   will_be_ice = F(Y+dY*dt(k),X+dX*dt(k))==1;

   tmp = land_backfill(:,:,k); 
   tmp(sometimes_ice & will_be_ice) = true; 

   %tmp = bwmorph(bwmorph(tmp,'close'),'open'); 
   if t(k)>=731747
      tmp(jakotrim & sometimes_ice) = false; 
   end
   tmp(always_land) = true; 
   tmp(never_land) = false;
   tmp = imfill(tmp,8,'holes'); 
   tmp(always_land) = true; 
   tmp(never_land) = false;

   land_backfill(:,:,k) = tmp; 
end

land_backfill = remove_icebergs(land_backfill,marine,always_land,never_land);

landsum(:,size(landsum,2)+1) = squeeze(sum(sum(land_backfill))); 

disp([datestr(now),': Backfill 1 finished'])

%% Backfill 2:24 months
% Repeat a similar back-filling process as above, but for 2 to 24 months. 
% Why is this code section separate from the one-month backfill section
% above? 
% 1. Beyond two or more consecutive months, we can let dt be 365/12 for
% every step, because there's not much variation in dt for multiple months. 
% 2. The step above would accurately account for moving ice if the velocity
% is at least one grid cell per month, whereas slow-moving ice may take a
% few months to advect an entire grid cell. 
% 3. We're doing our best to account for situations where a year of data
% might be missing. 

% Positions after 1 average month: 
X_adv = X + vtol*interp2(x,y,vx,X,Y)/12; 
Y_adv = Y + vtol*interp2(x,y,vy,X,Y)/12; 

for km = 2:24

    % Positions after km average months: 
    X_tmp = X_adv + vtol*interp2(x,y,vx,X_adv,Y_adv)/12; 
    Y_tmp = Y_adv + vtol*interp2(x,y,vy,X_adv,Y_adv)/12; 
    
    % We made tmp grids above bc we couldn't overwrite X_adv,Y_adv while solving for the new updated X_adv,Y_adv.   
    X_adv = X_tmp; 
    Y_adv = Y_tmp; 
    
    for k = (numel(t)-km):-1:1
       
       % Find out if it will be ice next month by interpolating next month's mask at today's grid point positions plus dt days of displacement: 
       %will_be_ice_old = interp2(x,y,single(land_backfill(:,:,k+km)),X_adv,Y_adv)==1;

       F.Values = single(flipud(land_backfill(:,:,k+km)));
       will_be_ice = F(Y_adv,X_adv)==1;
      
       tmp = land_backfill(:,:,k); 
       tmp(sometimes_ice & will_be_ice) = true; 
       if t(k)>=731747
          tmp(jakotrim & sometimes_ice) = false; 
       end
       tmp(always_land) = true; 
       tmp(never_land) = false;
       %tmp = bwmorph(bwmorph(tmp,'close'),'open'); 
       tmp = imfill(tmp,8,'holes'); 
       tmp(always_land) = true; 
       tmp(never_land) = false;

       land_backfill(:,:,k) = tmp; 
    end
    
    land_backfill = remove_icebergs(land_backfill,marine,always_land,never_land);

    disp([datestr(now),': Backfill ',num2str(km),' finished'])
end

clear X_adv Y_adv tmp k

landsum(:,size(landsum,2)+1) = squeeze(sum(sum(land_backfill))); 

% Write the data in case of crash: 
save(fn_backfill,'land_backfill','-v7.3')

%% Pre-final: Determine backfill vs obs
% Start with a backfill assumption, but overwrite wherever observations 
% are available. Work backwards, prioritizing the advection-corrected
% masking of next month's observations over the backfill cube we created
% above.

dt = diff(t); 

disp([datestr(now),': Starting Pre-final land cube.'])
   
for k = (length(t)-1):-1:1

   tmp = land_backfill(:,:,k); 

   % Find out if it will be ice next month by interpolating next month's mask at today's grid point positions plus dt days of displacement: 
   will_be_ice = interp2(x,y,single(land(:,:,k+1)),X+dX*dt(k),Y+dY*dt(k))==1;
   
   tmp(will_be_ice & sometimes_ice) = true; % If these pixels will move to a place that's ice next month, make them true today.
   
   % This month's terminus datapoints: 
   indp = T.t>=(t(k)-31) & T.t<=t(k); % Terminus points within the past 31 days
   indf = T.t<=(t(k)+31) & T.t>=t(k); % Terminus points within the next 31 days

   if any(indp) | any(indf)
      
      % Get terminus positions taken within 30 days in the past:  
      cxp = T.x(indp);
      cyp = T.y(indp);
      cpp = T.p(indp); 
      dtp = T.t(indp) - t(k); 
      
      % Account for <=30 days of advection:
      dxp = interp2(x,y,vx,cxp,cyp) .* dtp/365.25;
      dyp = interp2(x,y,vy,cxp,cyp) .* dtp/365.25;
      cxp = cxp - dxp; 
      cyp = cyp - dyp; 
      
      % Get terminus positions taken less than 30 days in the future:  
      cxf = T.x(indf);
      cyf = T.y(indf);
      cpf = T.p(indf); 
      dtf = T.t(indf) - t(k); 
      
      % Account for <=30 days of advection:
      dxf = interp2(x,y,vx,cxf,cyf) .* dtf/365.25;
      dyf = interp2(x,y,vy,cxf,cyf) .* dtf/365.25;
      cxf = cxf - dxf; 
      cyf = cyf - dyf; 
      
      % Loop through all priorities of data, starting with the lousiest: 
      for kp = 1:max(T.p)
         
         % Indices of advected past points that fall within unconstrained grid cells:
         indkpf = cpf==kp & cxf>x(1) & cxf<x(end) & cyf<y(1) & cyf>y(end); 



         % Upstream from future positions advected to today: 
         if any(indkpf)
            [xus,yus] = stream2_fast(x,yf,-vxf,-vyf,cxf(indkpf),cyf(indkpf),[step Npts]);
         else 
            xus = []; 
         end
         if ~isempty(xus)
            isf = isfinite(xus); 
            Us = gridbin(xus(isf),yus(isf),true(size(xus(isf))),x,y,@any); % gridded cells upstream of the terminus 
            Us = imfill(bwmorph(Us,'close'),4,'holes'); 
         else
            Us = false(size(X)); 
         end
         
         % Indices of advected past points that fall within unconstrained grid cells:
         indkpp = cpp==kp & cxp>x(1) & cxp<x(end) & cyp<y(1) & cyp>y(end); 
         


           % Downstream from past positions advected to today: 
           if any(indkpp)
               [xds,yds] = stream2_fast(x,yf,vxf,vyf,cxp(indkpp),cyp(indkpp),[step Npts]);
           else 
              xds = []; 
           end
           if ~isempty(xds)
               isf = isfinite(xds); 
               Ds = gridbin(xds(isf),yds(isf),true(size(xds(isf))),x,y,@any); % Grid cells that are downstream of the terminus 
               Ds = imfill(bwmorph(Ds,'close'),4,'holes'); 
           else
               Ds = false(size(X)); 
           end
         
         FillRegion = Us & ~rock; 
         CarveRegion = Ds & ~rock; 
         
         tmp(CarveRegion) = false; 
         tmp(FillRegion) = true; % overwrites any potential CarveRegion
   
         if t(k)>=731747
            tmp(jakotrim & sometimes_ice) = false; 
         end
         tmp(rock) = true; 
         tmp = imfill(tmp,8,'holes'); 
      
      end
   end
   
   land(:,:,k) = tmp; 
   disp([datestr(now),': Pre-final cube finished month ',num2str(k),' of ',num2str(length(t)),'.'])

end


land = remove_icebergs(land,marine,always_land,never_land);
landsum(:,size(landsum,2)+1) = squeeze(sum(sum(land_backfill))); 

%land_prefinal = land; 

% Write the data in case of crash: 
save(fn_prefinal,'land','-v7.3')

clear tmp indp indf Ds Us FillRegion CarveRegion st us ds dxf dyf cxd cyf cpf cpp cxp cyp dxp dyp 
%% Final: Determine FrontCarve vs obs
% Similar theory as above, but now we're moving forward. Start by assuming
% that if this month's grid call came from a place that was ocean last
% month, it can't suddenly become ice this month. Overwrite that assumption
% where observations are available. 

dt = diff(t); 

step = 0.2; % Number of streamline steps per grid cell.  
reach = 50e3; % meters of reach in each direction from measured termini
Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 

disp([datestr(now),': Starting Final land cube.'])
   
for k = 2:length(t)

   tmp = land(:,:,k); 

   % Find out if it will be ice next month by interpolating next month's mask at today's grid point positions plus dt days of displacement: 
   %was_ocean = interp2(x,y,single(land(:,:,k-1)),X-dX*dt(k-1),Y-dY*dt(k-1))==0;
   F.Values = single(flipud(land(:,:,k-1)));
   was_ocean = F(Y-dY*dt(k-1),X-dX*dt(k-1))==0;

   tmp(was_ocean & sometimes_ice) = false; % If these pixels came from a place that was definitely ocean a month ago, then they must still be ocean.  
   
   % This month's terminus datapoints: 
   indp = T.t>=(t(k)-31) & T.t<=t(k); % Terminus points within the past 31 days
   indf = T.t<=(t(k)+31) & T.t>=t(k); % Terminus points within the next 31 days
   
   % Get terminus positions taken within 30 days in the past:  
   cxp = T.x(indp);
   cyp = T.y(indp);
   cpp = T.p(indp); 
   dtp = T.t(indp) - t(k); 
   
   % Account for <=31 days of advection:
   dxp = interp2(x,y,vx,cxp,cyp) .* dtp/365.25;
   dyp = interp2(x,y,vy,cxp,cyp) .* dtp/365.25;
   cxp = cxp - dxp; 
   cyp = cyp - dyp; 
   
   % Get terminus positions taken less than 30 days in the future:  
   cxf = T.x(indf);
   cyf = T.y(indf);
   cpf = T.p(indf); 
   dtf = T.t(indf) - t(k); 
   
   % Account for <=31 days of advection:
   dxf = interp2(x,y,vx,cxf,cyf) .* dtf/365.25;
   dyf = interp2(x,y,vy,cxf,cyf) .* dtf/365.25;
   cxf = cxf - dxf; 
   cyf = cyf - dyf; 
   
   % Loop through all priorities of data, starting with the lousiest: 
   for kp = 1:max(T.p)
      
      % Indices of advected past points that fall within unconstrained grid cells:
      indkpf = cpf==kp & cxf>x(1) & cxf<x(end) & cyf<y(1) & cyf>y(end); 
      

          % Upstream from future positions advected to today: 
          if any(indkpf)
             [xus,yus] = stream2_fast(x,yf,-vxf,-vyf,cxf(indkpf),cyf(indkpf),[step Npts]);
          else
             xus = []; 
          end
         if ~isempty(xus)
            isf = isfinite(xus); 
            Us = gridbin(xus(isf),yus(isf),true(size(xus(isf))),x,y,@any); % gridded cells upstream of the terminus 
            Us = imfill(bwmorph(Us,'close'),4,'holes'); 
         else
            Us = false(size(X)); 
         end
      
%       % Indices of advected past points that fall within unconstrained grid cells:
       indkpp = cpp==kp & cxp>x(1) & cxp<x(end) & cyp<y(1) & cyp>y(end); 
%       

           % Downstream from past positions advected to today: 
           if any(indkpp)
               [xds,yds] = stream2_fast(x,yf,vxf,vyf,cxp(indkpp),cyp(indkpp),[step Npts]);
           else 
              xds = []; 
           end
           if ~isempty(xds)
               isf = isfinite(xds); 
               Ds = gridbin(xds(isf),yds(isf),true(size(xds(isf))),x,y,@any); % Grid cells that are downstream of the terminus 
               Ds = imfill(bwmorph(Ds,'close'),4,'holes'); 
           else
               Ds = false(size(X)); 
           end
         
      
      FillRegion = Us & ~rock; 
      CarveRegion = Ds & ~rock; 
      
      tmp(CarveRegion) = false; 
      tmp(FillRegion) = true; % overwrites any potential CarveRegion

      if t(k)>=731747
         tmp(jakotrim & sometimes_ice) = false; 
      end
      tmp(rock) = true; 
      tmp = imfill(tmp,8,'holes'); 
      tmp(rock) = true; 
   
   end
   
   land(:,:,k) = tmp; 
   disp([datestr(now),': Final cube finished month ',num2str(k),' of ',num2str(length(t)),'.'])

end

land = remove_icebergs(land,marine,always_land,never_land);
landsum(:,size(landsum,2)+1) = squeeze(sum(sum(land))); 

save(fn_final,'-v7.3')
disp 'ready for clean'

%%

vtol = 1.25; % velocity tolerance. 1.25 allows 25% velocity variation. 

land_preclean = land; 

dt = diff(t); 

disp([datestr(now),': Starting clean-up cube.'])
   
%for k = 2:(length(t)-1) 
%for k = [2:(length(t)-1) (length(t)-2):-1:2]
for k = (length(t)-1):-1:2

   tmp = land(:,:,k); 

   % Find out if it will be ice next month by interpolating next month's mask at today's grid point positions plus dt days of displacement: 
   
   F.Values = single(flipud(land(:,:,k+1)));
   will_be_ice = F(Y+dY*dt(k)*vtol,X+dX*dt(k)*vtol)==1;
   
   F.Values = single(flipud(land(:,:,k-1)));
   was_ice = F(Y-dY*dt(k-1)*vtol,X-dX*dt(k-1)*vtol)==1;

   tmp(sometimes_ice & will_be_ice & was_ice) = true; 
   tmp(sometimes_ice & ~was_ice) = false; % because the ice had to come from somewhere.

    if t(k)>=731747
        tmp(jakotrim & sometimes_ice) = false; 
    end
        tmp(rock) = true; 
        tmp = imfill(tmp,8,'holes'); 
   
   land(:,:,k) = tmp; 
   disp([datestr(now),': Final clean-up cube finished month ',num2str(k),' of ',num2str(length(t)),'.'])

end

clear tmp 

land = remove_icebergs(land,marine,always_land,never_land);
landsum(:,size(landsum,2)+1) = squeeze(sum(sum(land))); 



%% Clean 

vtol = 1.25; % velocity tolerance. 1.25 allows 25% velocity variation. 

dt = diff(t); 

disp([datestr(now),': Starting clean-up cube.'])
   
for k = [2:(length(t)-1) (length(t)-2):-1:2]

   tmp = land(:,:,k); 

   % Find out if it will be ice next month by interpolating next month's mask at today's grid point positions plus dt days of displacement: 
   
   F.Values = single(flipud(land(:,:,k+1)));
   will_be_ice = F(Y+dY*dt(k)*vtol,X+dX*dt(k)*vtol)==1;
   
   F.Values = single(flipud(land(:,:,k-1)));
   was_ice = F(Y-dY*dt(k-1)*vtol,X-dX*dt(k-1)*vtol)==1;

   tmp(sometimes_ice & will_be_ice & was_ice) = true; 
   tmp(sometimes_ice & ~was_ice & ~will_be_ice) = false; 

    if t(k)>=731747
        tmp(jakotrim & sometimes_ice) = false; 
    end
        tmp(rock) = true; 
        tmp = imfill(tmp,8,'holes'); 
   
   land(:,:,k) = tmp; 
   disp([datestr(now),': Nearly final clean-up cube finished month ',num2str(k),' of ',num2str(length(t)),'.'])

end

clear tmp 

land = remove_icebergs(land,marine,always_land,never_land);
landsum(:,size(landsum,2)+1) = squeeze(sum(sum(land))); 

disp([datestr(now),': Starting really truly final clean-up cube.'])
   
for k = length(t):-1:2

   tmp = land(:,:,k); 

   
   F.Values = single(flipud(land(:,:,k-1)));
   was_ice = F(Y-dY*dt(k-1)*vtol,X-dX*dt(k-1)*vtol)==1;

   tmp(sometimes_ice & ~was_ice) = false; 

        tmp(rock) = true; 
        tmp = imfill(tmp,8,'holes'); 
   
   land(:,:,k) = tmp; 
   disp([datestr(now),': Truly final clean-up cube finished month ',num2str(k),' of ',num2str(length(t)),'.'])

end

clear tmp 

land = remove_icebergs(land,marine,always_land,never_land);

% Write the data in case of crash: 
save(fn_clean,'land','-v7.3')

clear land_backfill land_prefinal marine never_land sometimes_ice tmp was_ocean will_be_ice th vx* vy*

%% Plot time series 

rocksum = squeeze(sum(sum(rock))); 

A_km2 = diff(x(1:2)/1000)^2; % grid cell area (km2) 

if ~devon

    col = cmocean('dense',size(landsum,2)+1); 

    figure
    hold on
    for k=1:size(landsum,2)

        plot(datetime(t,'convertfrom','datenum'),(landsum(:,k)-rocksum)*A_km2,'color',col(k+1,:),'linewidth',1)
    end
    axis tight
    box off
end

%%

ice = land; 
for k = 1:length(t)
    ice(:,:,k) = land(:,:,k).*(~rock); 
end
disp([datestr(now),' Saving NetCDF now...'])

clear land 
%% Save data 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Greenland monthly ice masks.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date_created',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Institution','NASA Jet Propulsion Laboratory (JPL), California Institute of Technology');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Please cite this dataset! Check the GitHub page for citation info.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'GitHub','https://github.com/chadagreene/greenland-coastlines');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'MATLAB_script',readme);

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
x_var_id = netcdf.defVar(ncid,'x','NC_SHORT',x_id);
netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate grid cell center');
netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,x_var_id,'units',        'meter');

% Define y
y_id     = netcdf.defDim(ncid,'y',length(y));
y_var_id = netcdf.defVar(ncid,'y','NC_SHORT',y_id);
netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate grid cell center');
netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,y_var_id,'units',        'meter');

% Define time
time_id     = netcdf.defDim(ncid,'time',length(t));
time_var_id = netcdf.defVar(ncid,'time','NC_INT',time_id);
netcdf.putAtt(ncid,time_var_id,'long_name',    'time');
netcdf.putAtt(ncid,time_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,time_var_id,'units',        'days since 1900-1-1 0:0:0');

% Define ice
ice_var_id = netcdf.defVar(ncid,'ice','NC_BYTE',[x_id y_id time_id]);
netcdf.putAtt(ncid, ice_var_id,'long_name','Binary grid indicates presence of ice.');
netcdf.putAtt(ncid, ice_var_id,'grid_mapping', 'mapping');

% Define rock
rock_var_id = netcdf.defVar(ncid,'rock','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid, rock_var_id,'long_name','Binary grid indicates presence of rock.');
netcdf.putAtt(ncid, rock_var_id,'grid_mapping', 'mapping');

% Compress and stop variable definition
netcdf.defVarDeflate(ncid,ice_var_id,true,true,9);
netcdf.defVarDeflate(ncid,rock_var_id,true,true,9);
netcdf.endDef(ncid);

%3. Place data
netcdf.putVar(ncid,x_var_id,x);
netcdf.putVar(ncid,y_var_id,y);
netcdf.putVar(ncid,time_var_id,t-datenum(1900,1,1,0,0,0));
netcdf.putVar(ncid,ice_var_id,uint8(ipermute(ice,[2 1 3])));
netcdf.putVar(ncid,rock_var_id,ipermute(uint8(rock),[2 1]));

%4. Close file 
netcdf.close(ncid)

disp([datestr(now),' All done'])


% This script analyzes Greenland ice masks that were created by
% greenland_termini_to_icemasks.m. 
% 
% Analysis includes 
% * overall secular change in flowline length from 1985 to 2021
% * overall secular change in ice sheet mass from 1985 to 2021 
% * seasonal length variability assessment 
% 
% Chad Greene, NASA/JPL, August 2022. 

%% Load data

filename = 'greenland_extruded_velocity_and_thickness.nc'; 

xf = double(ncread(filename,'x')); 
yf = double(ncread(filename,'y')); 

vx = double(permute(ncread(filename,'vx'),[2 1])); 
vy = double(permute(ncread(filename,'vy'),[2 1])); 
rock = permute(ncread(filename,'v_source'),[2 1])==0; 
th = double(permute(ncread(filename,'thickness'),[2 1])); 
ths = permute(ncread(filename,'thickness_source'),[2 1]);

[X,Y] = meshgrid(xf,yf); 

% Change thickness to hydrostatic where extrapolation? 
bed = bedmachine_interp('bed',X,Y,'greenland'); 
th2 = base2freeboard(bed)*9.3;
th3 = th; 
th3(ths==2 & th2>th) = th2(ths==2 &th2>th); 
clear th2 bed ths; 

% Ice masks: 
fn = '/Users/cgreene/Documents/GitHub/greenland-coastlines/data/greenland_monthly_ice_masks_2022-08-23.nc';
t = ncdateread(fn,'time'); 
ice = logical(permute(ncread(fn,'ice'),[2 1 3]));
x = double(ncread(fn,'x'));
y = double(ncread(fn,'y'));

%% Trim
% Trim the full ice sheet data to the region of the monthly mask file 
% (This is mainly for testing a small region). 

cols = ismember(xf,x); 
rows = ismember(yf,y); 
clear xf yf 

vx = vx(rows,cols); 
vy = vy(rows,cols); 
th = th(rows,cols); 
th3 = th3(rows,cols); 
rock = rock(rows,cols); 
X = X(rows,cols); 
Y = Y(rows,cols); 

clear rows cols keep 

%% Area and mass time series

A = abs(diff(x(1:2))*diff(y(1:2)));

M = 917*A*th3*1e-12; 

M_ts = squeeze(sum(sum(M.*double(ice)))); 

figure
plot(t,M_ts - M_ts(1))
box off 
axis tight
ylabel 'Cumulative mass change (Gt)'

%%

% Outlet glacier pixels: 
outlets = bwperim(imfill(all(ice,3) | rock,8,'holes'),8) & ~rock & ~bwperim(true(size(X)));  
%outlets = bwperim(imfill(all(ice,3) | rock,8,'holes'),8) & ~rock; % only do this to fill in edges (much slower)  
landterm =  bwperim(imfill(all(ice,3),8,'holes'),8) & ~outlets; 

xo = X(outlets); 
yo = Y(outlets); 

ds = stream2(x,y,vx,vy,xo,yo); 
dl = zeros(size(ds)); 

ind1985 = find(t==datetime(1985,3,15)); 
ind2021 = find(t==datetime(2021,3,15)); 

for k = 1:length(ds)
   xi = ds{k}(:,1); 
   yi = ds{k}(:,2); 
   d = pathdistpsn(xi,yi); 

   ice1 = interp2(x,y,ice(:,:,ind1985),xi,yi,'nearest');
   ice2 = interp2(x,y,ice(:,:,ind2021),xi,yi,'nearest');

   if all(ice1)
      f1 = numel(xi)+1; 
   else
      f1 = find(~ice1,1,'first'); 
   end
   l1 = d(f1-1);

   if all(ice2)
      f2 = numel(xi)+1; 
   else
      f2 = find(~ice2,1,'first'); 
   end
   l2 = d(f2-1);

   dl(k) = l2 - l1; 
end

figure
scatter(xo,yo,30,dl,'filled')
bedmachine('gl','greenland'); 
caxis([-1 1]*10000)
cmocean bal

%%

% Loop through each seed location to extrapolate terminus thickness along flowlines:  
V = ds; 
for k = 1:length(V) 
   V{k}(:,1) = dl(k); 
end

% Concatenate the cells: 
M = cell2mat(ds(:)); % 
V = cell2mat(V(:)); 

% Grid up the streamline data:  
isf = isfinite(V(:,1)) & isfinite(M(:,1)); 
Ds = gridbin(M(isf,1),M(isf,2),V(isf,1),x,y); 
holes = imfill(isfinite(Ds),8,'holes') & ~isfinite(Ds);

Ds = regionfill(Ds,holes); 

%%

us = stream2(x,y,-vx,-vy,xo,yo); 

% Loop through each seed location to extrapolate terminus thickness along flowlines:  
V = us; 
for k = 1:length(V) 
   V{k}(:,1) = dl(k); 
end

% Concatenate the cells: 
M = cell2mat(us(:)); % 
V = cell2mat(V(:)); 

% Grid up the streamline data:  
isf = isfinite(V(:,1)) & isfinite(M(:,1)); 
Dlg = gridbin(M(isf,1),M(isf,2),V(isf,1),x,y); 

Dlg(isnan(Dlg) & landterm) = 0; 

Dlg(isnan(Dlg) & any(ice,3)) = Ds(isnan(Dlg) & any(ice,3)); 

holes = imfill(isfinite(Dlg),8,'holes') & ~isfinite(Dlg) & ~rock;

Dlg = regionfill(Dlg,holes); 

figure
hdl = imagescn(x,y,Dlg/1000); 
%hdl.AlphaData = hdl.AlphaData*.8; 
axis image off
cmocean -bal
caxis([-1 1]*15)
cb = colorbar; 
ylabel(cb,'Glacier length change, 1985-2021 (km)')

[Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xlim,ylim,1000);
hold on
hh=image(xx,yy,Itmp(:,:,1:3));
uistack(hh,'bottom'); 
scalebarpsn('location','se')

%% Seasonal 

M_ts_lp = movmean(M_ts,12); 

[yr,mo,~] = datevec(t); 

keep = yr>=2000 & yr<=2018; 
yr = yr(keep); 
t = t(keep); 
ice = ice(:,:,keep); 
M_ts = M_ts(keep); 
M_ts_lp = M_ts_lp(keep); 
M_ts_hp = M_ts - M_ts_lp; 

%%
figure
plot(t,M_ts,'color',hex2rgb('#006DAA'))
hold on 
plot(t,M_ts_lp,'color',hex2rgb('#061A40'),'linewidth',1)
box off
axis tight

%%

uyr = unique(yr); 

col = parula(length(uyr)); 

figure
hold on
for k = 1:length(uyr)
   plot(datenum(t) - datenum(uyr(k),0,0),M_ts_hp,'color',col(k,:))
end
axis tight
ylabel 'Seasonal mass anomaly (Gt)'
Ms = climatology(M_ts_hp,t); 
plot([doy(t(1:12));doy(t(1:12))+365;doy(t(1:12))+365*2],repmat(Ms,3,1),'r','linewidth',2)

ft = sinefit(t,M_ts_hp); 
plot(1:365*2,sineval(ft,1:365*2),'linewidth',2,'color',rgb('green'))

set(gca,'xtick',[doy(t(1:12));doy(t(1:12))+365])
datetick('x','mmm','keepticks')
xlim([365 365*2])
ylim([-1 1]*18)
textcolorbar(uyr,'location','eo')

%%

x_ex = -180667;
y_ex = -2277772;

k = find(xo==x_ex & yo==y_ex)

xi = ds{k}(:,1); 
yi = ds{k}(:,2); 

d = pathdistpsn(xi,yi); 

lts = nan(numel(t),1); 
for kk=1:length(t)
   ind = interp2(x,y,ice(:,:,kk),xi,yi,'nearest');
   if all(ind)
      f = numel(xi)+1; 
   else
      f = find(~ind,1,'first'); 
   end
   lts(kk) = d(f-1);
end

lts_lp = movmean(lts,12); 

ti = t(1):t(end); 
lts_lpi = interp1(t,lts_lp,ti); 

ft = sinefit(t,lts-lts_lp); 

figure('pos',[20 20 300 200])
pl(1)=plot(t,lts/1000,'color',hex2rgb('#2ec4db'));
hold on
pl(2)=plot(t,lts_lp/1000,'color',hex2rgb('#3e0938'),'linewidth',.75);
pl(3)=plot(ti,(lts_lpi+sineval(ft,ti))/1000,'color',hex2rgb('#ec381d'));
box off
axis tight
set(gca,'fontsize',7)
ylabel('Flowline length anomaly (km)','fontsize',7)
txt=textcolorbar({'monthly mask';'12-mo mean';'12-mo mean + sinefit'},'colormap',[pl(1).Color;pl(2).Color;pl(3).Color],'fontsize',7,'location','ne');
xlim([datetime(2000,0,0) datetime(2018,12,31)])
txt.HorizontalAlignment = 'left'; 
txt.Position = [.68 1 0];
set(gca,'Ticklength',[.003 0.01])

yl = ylim; 
ylim([-.5 yl(2)])

txt2 = ntitle(['Amplitude: ',num2str(round(ft(1))),' m. Day of maximum length: ',datestr(ft(2),'mmmm dd'),'.'],'fontsize',6,'location','sw','color',pl(3).Color);

% export_fig('/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/mask_sinefit_example.jpg','-pdf','-r600','-p0.01','-painters')

%%

tmp = [];

A = zeros(size(ds)); 
ph = nan(size(ds)); 

parfor k = 1:length(ds)

   xi = ds{k}(:,1); 
   yi = ds{k}(:,2); 
   
   if numel(xi)>1
      d = pathdistpsn(xi,yi); 
   
      lts = nan(numel(t),1); 
      for kk=1:length(t)
         ind = interp2(x,y,ice(:,:,kk),xi,yi,'nearest');
         if all(ind)
            f = numel(xi)+1; 
         else
            f = find(~ind,1,'first'); 
         end
         lts(kk) = d(f-1);
      end

      lts_lp = movmean(lts,12); 

      ft = sinefit(t,lts-lts_lp); 
      A(k) = ft(1); 
      ph(k) = ft(2); 

   end

end

%%

[xtmp,ytmp] = pol2cart(ph*2*pi/365.25,A); 


% Loop through each seed location to extrapolate terminus thickness along flowlines:  
V = us; 
for k = 1:length(V) 
   V{k}(:,1) = xtmp(k); 
   V{k}(:,2) = ytmp(k); 
end

% Concatenate the cells: 
M = cell2mat(us(:)); % 
V = cell2mat(V(:)); 

% Grid up the streamline data:  
isf = isfinite(V(:,1)) & isfinite(V(:,2)) & isfinite(M(:,1)); 
Xtmp = gridbin(M(isf,1),M(isf,2),V(isf,1),x,y); 
Ytmp = gridbin(M(isf,1),M(isf,2),V(isf,2),x,y); 

holes = imfill(isfinite(Xtmp),8,'holes') & ~isfinite(Xtmp) & ~rock;

Xtmp = regionfill(Xtmp,holes); 
Ytmp = regionfill(Ytmp,holes); 

[Ph,AA] = cart2pol(Xtmp,Ytmp); 
Ph(Ph<0) = Ph(Ph<0)+2*pi; 
Ph = Ph*365.25/(2*pi); 

figure
hph = imagescn(x,y,Ph);
axis image off
hold on
caxis
caxis([0 365])
cmocean phase

alpha = AA/500; 
alpha(isnan(alpha)) = 0; 
alpha(alpha>1) = 1; 
hph.AlphaData = alpha; 


[Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xlim,ylim,1000);
hold on
hh=image(xx,yy,Itmp(:,:,1:3));
uistack(hh,'bottom'); 
set(gca,'pos',[0 0 1 1])

scalebarpsn('location','se')

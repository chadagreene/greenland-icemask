% This script loads Greenland terminus position shapefile data from Joughin's  
% annual winter terminus positions, and Goliber's TermPicks dataset (with 
% Cheng's CALFIN). This script also densifies each terminus line segment to
% ~24 m resolution. 
% 
% This entire script takes about 2 minutes to run on my laptop from 2019. 
% 
% Chad Greene, July 2022. 

grid_resolution = 120; % m final grid resolution (only used to determine grid_resolution/5 data density along terminus lines) 

%% Load Taryn Black's data 

B = shaperead('/Users/cgreene/Documents/data/coastlines/black/glacier_termini_v01.0.shp');
B = B([B.Quality_Fl]==0);

tb = datenum(datetime({B.SourceDate})); 
procOrderb = 7*ones(size(tb)); 
data_originb = 0*ones(size(procOrderb')); 

%% Load Measures terminus data (Joughin et al) 
% Note from Taryn Black: "I would rank flags 0 and 2 together as best, and
% flags 1 and 3 together as second-best. Flag 0 indicates that there were no
% issues with digitizing the terminus (i.e. I was confident about the position).
% Flag 2 rarely got used but essentially says "I wasn't sure about the terminus 
% so I checked another image and now I'm sure". Flag 1 indicates that I was 
% uncertain about the terminus position, usually due to an abundance of melange 
% at the terminus or illumination issues. Flag 3 is only used for Landsat-7 images
% when the SLC is off, and in those cases I digitized what I could see and 
% then connected straight across the missing strips in the images."

% Load Joughin dataset: 
tmp = load('/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/measures_greenland_terminus_picks_v2.mat'); 

% Trim Joughin dataset to only 0 and 2 flags:
J = tmp.S(ismember([tmp.S.Quality_Fl],[0 2])); 

% Trim away any Joughin data that's already in the Black dataset: 
redundant = ismember({J.Image_ID},{B.Image_ID});
J = J(~redundant); 
clear redundant

% Make a time array: 
tj = datenum(datetime({J.SourceDate})); 

procOrderj = 7*ones(size(tj)); % processing order, described below

clear tmp 

data_originj = 1*ones(size(procOrderj')); 

%% Load Goliber coastline data
% quality flags (and sums of each) 
% '0' or '00' (20639) no issues
% '01' (8) uncertainty due to clouds or other image issues 
% '02' (0) supplemented trace (two or more images) 
% '3' or '03' (1488) landsat 7 slc off
% '04' (0) incomplete/box method 
% '05' (15700) automatically assigned 
% '01, 04' (338) 
% '02, 04' (86) 
% '03, 04' (344) 
% '04, 04' (1) 
% '05, 04' (2315) 
% '10' (15522)
% '13' (4313)

fn = '/Users/cgreene/Documents/data/coastlines/TermPicks+CALFIN_V2/TermPicks+CALFIN_V2.shp';
G = shaperead(fn); % golibert dataset

% Get rid of redundant data (anything in Goliber's dataset that's also in Joughin's) 
G(ismember({G.ImageID},{J.Image_ID}) & strcmpi({G.Author},'black_taryn')) = []; 
tg = datenum(datetime([G.Year],[G.Month],[G.Day]));

% Set processing order (best data is highest number) 
procOrderg = zeros(size(tg),'uint8'); 
procOrderg(ismember({G(:).QualFlag}',{'04';'01, 04';'02, 04';'03, 04';'04, 04';'05, 04'})) = 1; % box method (likely least accurate) 
procOrderg(ismember({G(:).QualFlag}',{'0';'00'})) = 2; % manual, no issues 
procOrderg(strcmp({G(:).QualFlag}','05')) = 4; % automatically assigned gl id 
procOrderg(ismember({G(:).QualFlag}',{'3';'03';'13';'01';'02'})) = 3; % slc off or otherwise uncertain 
procOrderg(strcmp({G(:).QualFlag}','10')) = 6; % calfin, no issues

data_origing = 2*ones(size(procOrderg')); % TermPicks
data_origing(ismember({G(:).QualFlag},{'10';'13'})) = 3; % calfin

%% AutoTerm

load('/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/autoterm_reformatted_trimmed.mat','A','ta','sensor','err')

%keep = ismember(sensor,{'Landsat8','Sentinel1','Sentinel2'}) & err<100; 
keep = true(size(ta));

A = A(keep); 
ta = ta(keep); 

data_origina = 4*ones(size(ta)); 
procOrdera = 5*ones(size(ta)); 

clear sensor err keep
%% Concatenate terminus datasets 

% Remove unnecessary fields (to make them easily concatenatable): 
J = rmfield(J,{'Geometry','BoundingBox','Image_ID','Sensor','Glacier_ID','SourceDate','Quality_Fl'});
B = rmfield(B,{'Geometry','BoundingBox','Image_ID','Image_Tile','Glacier_ID','SourceDate','Quality_Fl'});
G = rmfield(G,{'Geometry','BoundingBox','GlacierID','Date','Year','Month','Day','DecDate','QualFlag','Satellite','ImageID','Author','Center_X','Center_Y'});

% Concatenate the datasets: 
C = cat(1,G,J,A,B); 
t = cat(2,tg,tj,ta',tb)';
procOrder = cat(2,procOrderg,procOrderj,procOrdera',procOrderb)';

data_origin = [data_origing;data_originj;data_origina;data_originb]; 

clear tg tj ta G J A fn procOrderg procOrderj procOrdera data_origing data_originj data_origina B procOrderb tb data_originb

%% When can we analyze changes most confidently? 
% Key takeaway: June or July 1985 to June 2021.

[yr,mo,~] = datevec(t); 

tm = datetime(1975,9,15):calmonths(1):datetime(2022,2,15); 
[yrtm,motm,~] = datevec(tm);
sm = nan(size(tm)); 
for k=1:length(tm)
    sm(k) = sum(yr==yrtm(k) & mo==motm(k)); 
end

figure
bar(tm,sm,'edgecolor','none')
box off
axis tight

clear yr mo tm yrtm motm sm 

%% Plot terminus position statistics 

[yr,~,~] = datevec(t); 
dy = doy(t); 

[sum(yr>=1972) sum(yr>=1985)]
% AutoTerm, CALFIN, TermPicks, MEaSURES v2:
[sum(yr>=1972 & data_origin==4) sum(yr>=1972 & data_origin==3) sum(yr>=1972 & data_origin==2) sum(yr>=1972 & data_origin==1)]
[sum(yr>=1985 & data_origin==4) sum(yr>=1985 & data_origin==3) sum(yr>=1985 & data_origin==2) sum(yr>=1985 & data_origin==1)]


yr_edge = min(yr)-.5:max(yr)+.5;
dy_edge = .5:365.5;

for k = 0:4
   counts_yr(k+1,:) = histcounts(yr(data_origin==k),yr_edge);
   counts_dy(k+1,:) = histcounts(dy(data_origin==k),dy_edge); 
end

% ytickcolor = hex2rgb('#b6b39d'); 
% ytickcolor = hex2rgb('#b1ae95'); 
% xcolor = hex2rgb('#093824'); 
% hc = hex2rgb('#BF4E30'); 
% hc = hex2rgb({'#406E8E';'#522b47';'#BF4E30';'#F4D8CD'}); 

ytickcolor = hex2rgb('#9A9B99'); 
xcolor = hex2rgb('#3E3F3D'); 
hc = hex2rgb({'#AA100E';'#F25F5C';'#247BA0';'#FFE066';'#70C1B3'}); 

fs = 7; % fontsize 

figure('position',[100 100 320 280])
ax(1)=subsubplot(2,1,1,'vpad',0.1);
hb = bar(min(yr):max(yr),counts_yr,1,'stacked','edgecolor','none');
hb(1).FaceColor = hc(1,:); 
hb(2).FaceColor = hc(2,:); 
hb(3).FaceColor = hc(3,:); 
hb(4).FaceColor = hc(4,:); 
hb(5).FaceColor = hc(5,:); 
axis tight 
box off
set(gca,'xcolor',xcolor,'fontsize',fs)

h1 = hline([2000 5000:5000:25000]); 
set(h1(:),'color',ytickcolor,'linewidth',0.25); 
uistack(h1,'bottom')
xlim([1970 2022])
xl = xlim; 
txt1 = text(xl(1)*ones(size([2000 5000:5000:20000])),[2000 5000:5000:20000],num2str(([2 5:5:20])'),...
   'fontsize',fs,'color',ytickcolor,'horiz','left','vert','bot');
txt1(7) = text(xl(1),25000,'25 thousand terminus positions per year',...
   'fontsize',fs,'color',ytickcolor,'horiz','left','vert','bot');
set(gca,'ycolor','none','xcolor',xcolor); 

ax(2) = subsubplot(2,1,2,0.1); 
hb = bar(1:365,counts_dy,1,'stacked','edgecolor','none');
hb(1).FaceColor = hc(1,:); 
hb(2).FaceColor = hc(2,:); 
hb(3).FaceColor = hc(3,:); 
hb(4).FaceColor = hc(4,:); 
hb(5).FaceColor = hc(5,:); 
axis tight 
box off
set(gca,'xtick',datenum(0,1:12,15),'color','none','fontsize',fs,'xcolor',xcolor)
datetick('x','keepticks')
xlim([0 366])

h2 = hline(200:200:1200); 
set(h2(:),'color',ytickcolor,'linewidth',0.25); 
uistack(h2,'bottom')
xl = xlim; 
txt2 = text(xl(1)*ones(size(200:200:1000)),200:200:1000,num2str((200:200:1000)'),...
   'fontsize',fs,'color',ytickcolor,'horiz','left','vert','bot');
txt2(6) = text(xl(1),1200,'1200 total terminus positions per day',...
   'fontsize',fs,'color',ytickcolor,'horiz','left','vert','bot');
set(gca,'ycolor','none','xcolor',xcolor); 

ht1 = textcolorbar({'Zhang';'Cheng';'Goliber';'Joughin v2';'Black'},...
   'colormap',flipud(hc),'location','ne','fontsize',fs,'background','w'); 
ht1.Margin = 0.01; 
ht1.Position(2)=.94;

set(gcf,'renderer','painters')
% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/greenland-icemask/figures/terminus_histograms.jpg','resolution',600)

% housekeeping
clear yr mo dy ax h1 h2  histcolor  fs xl txt*

%% Densify the terminus segments: 

animateDensification = false; 
saveExamplePlot = false; 

for k=1:length(C)
   tmpx = [C(k).X]; 
   tmpy = [C(k).Y]; 

   [tmpxc,tmpyc] = polysplit(tmpx,tmpy);

   for kc = 1:length(tmpxc)
      % add small random numbers (<1 m variation) to ensure all datapoints are unique
      if numel(tmpxc{kc})>1
         [tmpxc{kc},tmpyc{kc}] = psnpath(tmpxc{kc}+0.1*rand(size(tmpxc{kc})),tmpyc{kc},grid_resolution/5); % densify to a resolution of 5 points per target grid cell. 
      end
   end

   [tmpx,tmpy] = polyjoin(tmpxc,tmpyc); 
   
   % This creates an animation to check densification as it goes: 
   if animateDensification
      clf
      scatter([C(k).X],[C(k).Y],40,1:length([C(k).X]),'filled')
      axis image
      hold on
      plot(tmpx,tmpy,'.')
      xlabel 'easting (m)'
      ylabel 'northing (m)' 
      cb = colorbar; 
      ylabel(cb,'datapoint index')
      drawnow
      pause(.1)
   end
   
   if (k==5 & saveExamplePlot)
      % 263 points of k=5 follow discontinuous path with two line segments.
      figure('position',[100 100 362 271])
      scatter([C(k).X]/1000,[C(k).Y]/1000,10,1:length([C(k).X]),'filled')
      axis image
      hold on
      plot(tmpx/1000,tmpy/1000,'.','markersize',0.5)
      xlabel('easting (km)','fontsize',7)
      ylabel('northing (km)','fontsize',7)
      cb = colorbar; 
      ylabel(cb,'datapoint index','fontsize',7)
      set(gca,'fontsize',7)
      drawnow
      export_fig('/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/terminus_densification.jpg','-painters','-r600','-pdf')
      
   end
   
   C(k).X = tmpx'; 
   C(k).Y = tmpy'; 
   C(k).t = repmat(t(k),size(C(k).X)); 
   C(k).procOrder = repmat(procOrder(k),size(C(k).X)); 
   C(k).do = repmat(data_origin(k),size(C(k).X)); 

   if mod(k,10000)==0
      [k round(100*k/length(C))]
   end
end

x = single([C.X]);
y = single([C.Y]);
t = single([C.t]);
p = [C.procOrder];
do = [C.do]; 

% Trim to finite data: 
isf = isfinite(x); 
x = x(isf); 
y = y(isf); 
t = t(isf); 
p = p(isf); 
do = do(isf); 

% Trim slow-moving data: 
[mask,x_bm,y_bm] = bedmachine_data('mask','greenland');
bad = interp2(x_bm,y_bm,double(mask==1),x,y,'linear')==1; % linear interp of double mask ensures it's actually rock, not just nearest neighbor
x(bad) = [];
y(bad) = [];
t(bad) = [];
p(bad) = [];
do(bad) = []; 

% housekeeping: 
clear procOrder tmp* saveExamplePlot animateDensification k kc C isf bad

%% 

readme = 'Compiled TermPicks+CALFIN and Joughin annual terminus positions v2, with a culled version of AutoTerm, and Taryn Blacks recent 6-day and monthly positions, all densified to 24 m resolution. This data was compiled by terminus_data_densifier. t is datenum, x & y are north polar stereogrpahic (ps70) meters, and p is the priority assignment based on assumed data quality, with highest numbers being the best data.'; 

% save('/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/terminus_data_densified_2023-01-09.mat','x','y','t','p','readme')

%% Visual inspection 

[mask,x_bm,y_bm] = bedmachine_data('mask','greenland');

ind = 1:2:length(x); 
try
    yr = year(double(t)); 
catch
    [yr,~,~,] = datevec(double(t)); 
end

figure('pos',[ 3    57   888   890])
subplot(1,2,1)
set(gca,'pos',[0 0 .5 1])
maskoverlay(x_bm,y_bm,mask==1,'color',rgb('gray'));
hold on
axis image off
fastscatter(x(ind),y(ind),yr(ind),'markersize',1)
cmocean thermal
colorbar('location','south')
caxis([1985 2022])
ax = gca; 

subplot(1,2,2)
set(gca,'pos',[0.5 0 .5 1])
hh= maskoverlay(x_bm,y_bm,mask==1,'color',rgb('gray'));
hold on
axis image off
fastscatter(x(ind),y(ind),single(p(ind)),'markersize',1)
ax(2) = gca; 
label_greenland_glaciers('fontsize',9)

cb = colorbar('south');
set(cb,'xtick',1:7,'xticklabel',{'box';'manual';'slc off';'auto id';'AutoTerm';'CALFIN';'MEaSUREs'},'fontsize',7)
caxis([.5 7.5])
colormap(gca,parula(7))
linkaxes(ax,'xy')
return
%% Repeat above for an example image 

figa = true; % figure a or b 

yr = year(double(t)); 

figure('pos',[ 3  406 617 520])
subplot(1,2,1)
set(gca,'pos',[0 0 .5 1])
hold on
fastscatter(x(do==4),y(do==4),yr(do==4),'markersize',1)
fastscatter(x(do==2),y(do==2),yr(do==2),'markersize',1)
fastscatter(x(do==3),y(do==3),yr(do==3),'markersize',1)
fastscatter(x(do==1),y(do==1),yr(do==1),'markersize',1)
axis image off
cmocean thermal
cb = colorbar('location','north');
clim([1985 2022])
label_greenland_glaciers('fontsize',9)
ax = gca; 
if figa
set(cb,'fontsize',9,'XCOLOR','k') 
else
set(cb,'fontsize',9,'XCOLOR','w') 
end

subplot(1,2,2)
set(gca,'pos',[0.5 0 .5 1])
hold on
plot(x(do==4),y(do==4),'.','color',hc(4,:),'markersize',1)
plot(x(do==2),y(do==2),'.','color',hc(2,:),'markersize',1)
plot(x(do==3),y(do==3),'.','color',hc(3,:),'markersize',1)
plot(x(do==1),y(do==1),'.','color',hc(1,:),'markersize',1)
axis image off
ax(2) = gca; 
label_greenland_glaciers('fontsize',9)

cb = colorbar('north');
clim([.5 4.5])
set(cb,'xtick',1:4,'xticklabel',{'MEaSUREs v2';'TermPicks';'CALFIN';'AutoTerm'},'fontsize',9)
if figa
set(cb,'fontsize',9,'XCOLOR','k') 
else
set(cb,'fontsize',9,'XCOLOR','w') 
end
colormap(gca,hc)
linkaxes(ax,'xy')

if figa
% figure a: 
axis([ -268869     -228154    -2001612    -1928509]) 
else
% figure b: 
axis([84827      113953    -3045741    -2996092]) 
end

[Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xlim,ylim,1000);
h = image(xx,yy,Itmp(:,:,1:3));
uistack(h,'bottom'); 
axes(ax(1))
h = image(xx,yy,Itmp(:,:,1:3));
uistack(h,'bottom'); 

if figa
scalebarpsn('location','sw','color','w')
%exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/terminus_data_densifier_map_a.jpg','resolution',600)
else
scalebarpsn('location','sw')
%exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/terminus_data_densifier_map_b.jpg','resolution',600)

end

%%

T = load('terminus_data_densified_2023-01-09.mat','t'); 

spacing_m = 120/5; % along-path spacing of the densified terminus data is 5 points per grid cell. 
t = datenum(datetime(1972,9,15):calmonths(1):datetime(2022,2,15));
[yr,mo,~] = datevec(t); 

[yr_data,mo_data,~] = datevec(double(T.t)); 

obs_km = nan(size(t)); 
for k=1:length(t)
    obs_km(k) = sum(yr_data==yr(k) & mo_data==mo(k))*spacing_m/1000; 
end

figure
bar(datetime(t,'convertfrom','datenum'),obs_km)
axis tight
box off
ylabel('Terminus observation data per month (km)')


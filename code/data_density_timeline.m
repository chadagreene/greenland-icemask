

%% Load data 

load('icemask_catchment_analysis_2023-01-24.mat','names') 
names{244} = 'Ab Drachmann L Bistrup'; % just abbreviating for display purposes


% Ice masks: 
fn = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_monthly_ice_masks_2023-03-22.nc';
x = double(ncread(fn,'x'));
y = double(ncread(fn,'y'));
t = ncdateread(fn,'time');
%t = datetime(1972,9,15):calmonths(1):datetime(2022,02,15);

obs = ncread(fn,'observation_data')';

%% Plot data 
cm = cmocean('dense'); 
s1 = sum(obs); 
s2 = sum(obs'); 

fs = 7;
%dind = [1 34:60:590 594]; 
dind = [ 34:60:590 ]; 

txtcol = 0.4*[1 1 1];

figure

ax = axes;
ax.Position = [.05 .1 .6 .7];
imagescn(datenum(t),1:261,obs)
cmocean dense
set(gca,'colorscale','log')
caxis([1 500])
disableDefaultInteractivity(gca)
axis ij
set(gca,'ytick',0:10:300,'xtick',datenum(t(dind)),'xticklabel',datestr(t(dind),'yyyy'),'fontsize',fs)
ylabel('Catchment ID','fontsize',fs)
ylim([0.5 261])

ax(2) = axes('position',[.05 .8 .6 .15]); 
bar(t,s1/1000,1,'facecolor',cm(end,:),'edgecolor','none','facealpha',1)
box off
axis tight 
set(gca,'xtick',[],'fontsize',fs)
ylim([0 15])
ntitle({' Total monthly observations';' (thousands of km)'},'location','nw','fontsize',fs)


ax(3) = axes('position',[.65 .1 .15 .7]); 
barh(1:261,s2/1000,1,'facecolor',cm(end,:),'edgecolor','none','facealpha',1)
box off
axis tight 
set(gca,'ytick',[],'fontsize',fs,'xcolor','none')
xlim([0 40])
axis ij
ylim([0.5 261])

%ntitle(' Total monthly observations (km)','location','nw','fontsize',fs)

% txt(1) = text(s2(214)/1000,214,names{214},'horiz','left','vert','mid','fontangle','italic','color',txtcol,'fontsize',fs-1); 
% txt(2) = text(s2(218)/1000,218,names{218},'horiz','left','vert','mid','fontangle','italic','color',txtcol,'fontsize',fs-1); 

tmp = [6 38 50 63 86 102 139 190 196 214 218 246 249];
for k = 1:length(tmp)
    txt(k) = text(s2(tmp(k))/1000,tmp(k),names{tmp(k)},'horiz','left','vert','mid','fontangle','italic','color',txtcol,'fontsize',fs-1); 
end

export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/datadensity_timeline.jpg','-r900','-p0.01')

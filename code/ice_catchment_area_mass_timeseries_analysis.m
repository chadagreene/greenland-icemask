
load('icemask_catchment_analysis_1972-2022_v1.mat') 
load('greenland_ice_masks_1972-2022_v1_ice_sum')
load('testall.mat','lt') % land terminating 

A_err = A_err'; 
A_ts = A_ts'; 
M_err_pick_ts = M_err_pick_ts'; 
M_err_th_ts = M_err_th_ts'; 
M_ts = M_ts'; 

M_ts_err = hypot(M_err_th_ts,M_err_pick_ts); 


A_interann = movmean(A_ts,12,2); 
A_seasonal = A_ts - A_interann; 

M_interann = movmean(M_ts,12,2); 
M_seasonal = M_ts - M_interann; 

%t = datetime(t,'convertfrom','datenum'); 

A_tot = sum(A_ts); 
A_tot_seasonal = sum(A_seasonal); 
A_tot_interannual = sum(A_interann); 
A_tot_err = rssq(A_err); 

M_tot = sum(M_ts); 
M_tot_seasonal = sum(M_seasonal); 
M_tot_interannual = sum(M_interann); 
M_tot_err = rssq(M_ts_err); 

observation_data = ncread('greenland_ice_masks_1972-2022_v1.nc','observation_data');
x = ncread('greenland_ice_masks_1972-2022_v1.nc','x');
y = ncread('greenland_ice_masks_1972-2022_v1.nc','y');
catchment = permute(ncread('greenland_ice_masks_1972-2022_v1.nc','catchment'),[2 1]);

%%

[yr,mo,dy] = datevec(t); 

yrs = 2014:2020; 

Ms = nan(12,length(yrs)); 
As = Ms;

for k = 1:12
    ind = yr>=yrs(1) & yr<=yrs(end) & mo==k; 
    As(k,:) = A_tot_seasonal(ind); 
end

As_median = median(As,2); 
As_std = std(As,[],2); 

for k = 1:12
    ind = yr>=yrs(1) & yr<=yrs(end) & mo==k; 
    Ms(k,:) = M_tot_seasonal(ind); 
end

Ms_median = median(Ms,2); 
Ms_std = std(Ms,[],2); 

%%
% 
% for k = 1:12
%     As(k,isoutlier(As(k,:))) = nan; 
%     Ms(k,isoutlier(Ms(k,:))) = nan; 
% end
% 
% As_mean = mean(As,2,'omitnan'); 
% As_error = std(As,[],2,'omitnan');  
% 
% Ms_mean = mean(Ms,2,'omitnan'); 
% Ms_error = std(Ms,[],2,'omitnan');  
% 
% figure
% subplot(2,1,1) 
% boundedline(1:12,Ms_median,Ms_std)
% box off 
% axis tight
% 
% subplot(2,1,2) 
% boundedline(1:12,Ms_mean,Ms_error)
% box off 
% axis tight

%% 

date_range = datenum(t)>=datenum('jan 1, 2000'); 

pvA = polyfit(datenum(t(date_range)),A_tot(date_range),1); 
A_rate_km2_per_yr = pvA(1)*365.25/1e6

pvM = polyfit(datenum(t(date_range)),M_tot(date_range),1); 
M_rate_Gt_per_yr = pvM(1)*365.25

%% 

col1 = hex2rgb('#2A2B2A'); 
col2 = hex2rgb('#FF8966'); 
col3 = hex2rgb('#E5446D'); 


col1 = hex2rgb('#717EC3'); 
col2 = hex2rgb('#EE8434'); 
col3 = hex2rgb('#C95D63'); 

axcol = 0.1*[1 1 1]; 

fs = 7;


col = parula(length(yrs)+1); 
%col = cmasher('ember',length(yrs)+2); 
col = col(1:end-1,:); 

figure('pos',[139.00        474.00        323.00*.89*2*1.06        420.00*.8])

subsubplot(2,2,1,'vpad',0.02,'hpad',0.06)
[hb1,hb2]=boundedline(datenum(t),A_tot/1e6 -A_tot(end)/1e6,A_tot_err/1e6,'nan','gap','color',col1,'alpha','transparency',0.15);
hold on
p2=plot(datenum(t),A_tot_interannual/1e6-A_tot(end)/1e6,'color',col2);
tmp = polyval(pvA,datenum(t(date_range)))/1e6-A_tot(end)/1e6;
p3=plot(datenum(t(date_range)),tmp,'color',col3);
box off 
axis tight
xlim([datenum('jul 1, 1985') datenum('mar 1, 2022')])
datetick('x','yyyy','keeplimits')
set(gca,'fontsize',fs,'ycolor',axcol,'color','none')
ylabel('Ice sheet area (km^2) w.r.t. 2022','fontsize',fs) 
ylim(ylim + 100*[-1 1])
t1 = datenum('jan 15, 2014'); 
tmp = polyval(pvA,t1)/1e6-A_tot(end)/1e6;
text(t1,tmp,[num2str(round(A_rate_km2_per_yr)),' km^{2} yr^{-1}'],'horiz','left','vert','bot','color',col3,'fontsize',6)

leg = legend([hb1,p2,p3],{'Ice sheet monthly total','12 month moving mean','2000-2022 trend'},'location','southwest');
legend boxoff 
leg.ItemTokenSize(1)=5;
ntitle(' a ','fontsize',8,'fontweight','bold','location','nw')

subsubplot(2,2,3,'vpad',0.02,'hpad',0.06)
boundedline(datenum(t),M_tot-M_tot(end),M_tot_err,'nan','gap','color',col1,'alpha','transparency',0.15)
hold on
plot(datenum(t),M_tot_interannual-M_tot(end),'color',col2)
p3=plot(datenum(t(date_range)),polyval(pvM,datenum(t(date_range)))-M_tot(end),'color',col3);
box off 
axis tight
xlim([datenum('jul 1, 1985') datenum('mar 1, 2022')])
datetick('x','yyyy','keeplimits')
set(gca,'fontsize',fs,'ycolor',axcol,'color','none')
ylabel('Ice sheet mass (Gt) w.r.t. 2022','fontsize',fs) 
ylim(ylim + 30*[-1 1])
tmp = polyval(pvM,t1)-M_tot(end);
text(t1,tmp,[num2str(round(M_rate_Gt_per_yr)),' Gt yr^{-1}'],'horiz','left','vert','bot','color',col3,'fontsize',6)
ntitle(' b ','fontsize',8,'fontweight','bold','location','nw')


subsubplot(2,2,2,'vpad',0.02,'hpad',0.06)
hold on
[hb1,hb2] = boundedline(-11:24,repmat(As_median,3,1)/1e6,repmat(As_std/1e6,3,1),'color',col1,'alpha','transparency',0.15); 
for k = 1:length(yrs)
    ind = yr==yrs(k); 
    plot(1:12,A_tot_seasonal(ind)/1e6,'.','color',col(k,:),'markersize',6)
end

axis([0.5 12.5 -180 180])
set(gca,'xtick',1:12,'xticklabel',datestr(datenum(0,1:12,15),'mmm'),'fontsize',fs,'ycolor',axcol,'color','none')
ylabel('Seasonal area anomaly (km^2)','fontsize',fs) 
textcolorbar(yrs,'colormap',col,'location','ne','fontsize',fs)
ntitle(' c ','fontsize',8,'fontweight','bold','location','nw')


subsubplot(2,2,4,'vpad',0.02,'hpad',0.06)
hold on
[hb1,hb2] = boundedline(-11:24,repmat(Ms_median,3,1),repmat(Ms_std,3,1),'color',col1,'alpha','transparency',0.15); 

for k = 1:length(yrs)
    ind = yr==yrs(k); 
    plot(1:12,M_tot_seasonal(ind),'.','color',col(k,:),'markersize',6)
end

axis([0.5 12.5 -48 48])
set(gca,'xtick',1:12,'xticklabel',datestr(datenum(0,1:12,15),'mmm'),'fontsize',fs,'xcolor',axcol,'ycolor',axcol,'color','none')
ylabel('Seasonal mass anomaly (Gt)','fontsize',fs) 
xtickangle(0)
ntitle(' d ','fontsize',8,'fontweight','bold','location','nw')

% export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_seasonal_area_mass_anomalies.jpeg','-r900','-p0.01','-painters')
% export_fig('/Users/cgreene/Documents/papers/greene2023greenland/figures/jpg/fig_ED01.jpg','-r900','-p0.01')
% exportgraphics(gcf,'/Users/cgreene/Documents/papers/greene2023greenland/figures/eps/fig_ED01.eps','ContentType','vector')
return
%% Highlight slide 

figure
boundedline(datenum(t),M_tot-M_tot(end),M_tot_err,'nan','gap','color',col1,'alpha','transparency',0.15)
hold on
%plot(datenum(t),M_tot_interannual-M_tot(end),'color',col2)
p3=plot(datenum(t(date_range)),polyval(pvM,datenum(t(date_range)))-M_tot(end),'color',col3);
box off 
axis tight
xlim([datenum('jul 1, 1985') datenum('mar 1, 2022')])
datetick('x','yyyy','keeplimits')
set(gca,'fontsize',fs,'ycolor',axcol,'color','none')
ylabel('Ice sheet mass (Gt) w.r.t. 2022','fontsize',fs) 
ylim(ylim + 30*[-1 1])
tmp = polyval(pvM,t1)-M_tot(end);
text(t1,tmp,[num2str(round(M_rate_Gt_per_yr)),' Gt yr^{-1}'],'horiz','left','vert','bot','color',col3,'fontsize',6)

set(gcf,'pos',[560   612   378   236])

export_fig highlight_slide_2.jpg -r1200 -p0.01

%%

cmap = cmocean('thermal',5); 
col3 = cmap(3,:); 

figure
[hb1,hb2]=boundedline(datenum(t),M_tot-M_tot(end),M_tot_err,'nan','gap','color',col1,'alpha','transparency',0.15)
hold on
%plot(datenum(t),M_tot_interannual-M_tot(end),'color',col2)
p3=plot(datenum(t(date_range)),polyval(pvM,datenum(t(date_range)))-M_tot(end),'color',col3);
box off 
axis tight
xlim([datenum('jul 1, 1985') datenum('mar 1, 2022')])
datetick('x','yyyy','keeplimits')
set(gca,'fontsize',fs,'ycolor',axcol,'color','none')
ylabel('Ice lost to retreat (Gt)','fontsize',fs) 
ylim(ylim + 30*[-1 1])
tmp = polyval(pvM,t1)-M_tot(end);
text(t1,tmp,[num2str(round(M_rate_Gt_per_yr)),' Gt yr^{-1}'],'horiz','left','vert','bot','color',col3,'fontsize',8,'fontweight','bold')

set(gcf,'pos',[560   612   378   236],'color','k')
set(gca,'xcolor',.8*[1 1 1],'ycolor',.8*[1 1 1])

hb1.Color = cmap(5,:); 
hb2.FaceColor = cmap(5,:); 
hb2.FaceAlpha = 0.2;
hb1.LineWidth = 1; 

p3.LineWidth = 1; 
set(gca,'fontsize',8)
set(gca,'XTicklabelrotation',0)
set(gca,'ytick',0:100:1000)

% export_fig highlight_slide_2b.jpg -r1200 -p0.01

%%

m1 = 5; % month of seasonal max 
m2 = 9; % month of seasonal min

[diff(As_median([m1 m2]))/1e6 hypot(As_std(m1),As_std(m2))/1e6 diff(Ms_median([m1 m2])) hypot(Ms_std(m1),Ms_std(m2))]


%%

ind1 = 156; 
ind2 = 594; 

A_diff_km2 = (A_interann(:,ind2) - A_interann(:,ind1))/1e6;
A_diff_km2_err = hypot(A_err(:,ind1),A_err(:,ind2))/1e6;

round([sum(A_diff_km2) rssq(A_diff_km2_err)])

M_diff_Gt = (M_interann(:,ind2) - M_interann(:,ind1)); 
M_diff_Gt_error = hypot(M_ts_err(:,ind1),M_ts_err(:,ind2));

round([sum(M_diff_Gt) rssq(M_diff_Gt_error)])

[~,inda] = sort(A_diff_km2,'ascend');
[~,indm] = sort(M_diff_Gt,'ascend');

N = 5; 
names(indm(1:N))
round([A_diff_km2(indm(1:N)) A_diff_km2_err(indm(1:N))])
round([M_diff_Gt(indm(1:N)) M_diff_Gt_error(indm(1:N))])


names(inda(1:N))
round([A_diff_km2(inda(1:N)) A_diff_km2_err(inda(1:N))])
round([M_diff_Gt(inda(1:N)) M_diff_Gt_error(inda(1:N))])

ind = 21; 
names(ind)
([A_diff_km2(ind) A_diff_km2_err(ind)])
([M_diff_Gt(ind) M_diff_Gt_error(ind)])

figure('pos',[ 55.00        618.00        392.00*3.5/4.9        258.00*3.5/4.9])

axis([-1000 20 -220 10])
hold on
hline(0,'color',.8*[1 1 1])
vline(0,'color',.8*[1 1 1])
plot([A_diff_km2 A_diff_km2]',M_diff_Gt' + 0.5.*[-M_diff_Gt_error M_diff_Gt_error]','-','color',rgb('blue'))
hold on
plot(A_diff_km2' + .5.*[-A_diff_km2_err A_diff_km2_err]' ,[M_diff_Gt M_diff_Gt]','-','color',rgb('blue'))
plot(A_diff_km2,M_diff_Gt,'.','markersize',5,'color',rgb('blue'))
nd = abs(A_diff_km2)>50;
text(A_diff_km2(nd),M_diff_Gt(nd),names(nd),'color',.4*[1 1 1],'fontsize',5,'horiz','center','vert','bot')
xlabel('Area change 1985-2022 (km^2)','fontsize',7)
ylabel('Mass change 1985-2022 (Gt)','fontsize',7)
set(gca,'fontsize',7)

%% 
has_data  = sum(observation_data)>0;

sum(has_data' & abs(M_diff_Gt)>M_diff_Gt_error & abs(A_diff_km2)>A_diff_km2_err)

%% 

M_seasonal_median = nan(261,12); 
A_seasonal_median = nan(261,12); 
M_seasonal_mad = nan(261,12); 
A_seasonal_mad = nan(261,12); 

ind = 1:12:84;

for k = 1:12
   A_seasonal_median(:,k) = median(A_seasonal(:,495+ind+k),2);
   M_seasonal_median(:,k) = median(M_seasonal(:,495+ind+k),2);
   A_seasonal_mad(:,k) = mad(A_seasonal(:,495+ind+k),1,2);
   M_seasonal_mad(:,k) = mad(M_seasonal(:,495+ind+k),1,2);
end



[A_seasonal_max,A_max_month] = max(A_seasonal_median,[],2);
[M_seasonal_max,M_max_month] = max(M_seasonal_median,[],2);
[A_seasonal_min,A_min_month] = min(A_seasonal_median,[],2);
[M_seasonal_min,M_min_month] = min(M_seasonal_median,[],2);

A_seasonal_range = A_seasonal_max-A_seasonal_min;
M_seasonal_range = M_seasonal_max-M_seasonal_min;

for k = 1:261
    A_seasonal_error(k,:) = 1.4826*hypot(A_seasonal_mad(k,A_max_month(k)),A_seasonal_mad(k,A_min_month(k)));
    M_seasonal_error(k,:) = 1.4826*hypot(M_seasonal_mad(k,M_max_month(k)),M_seasonal_mad(k,M_min_month(k)));
end

seasonal_area = has_data' & A_seasonal_range>A_seasonal_error; 
seasonal_area_pct = sum(seasonal_area)*100/sum(has_data)

seasonal = has_data' & A_seasonal_range>1*A_seasonal_error & M_seasonal_range>1*M_seasonal_error & M_seasonal_range>1e-5; 

N=10; 
[~,ind] = sort(M_seasonal_range,'descend');
names(ind(1:N))
[A_seasonal_range(ind(1:N))/1e6 A_seasonal_error(ind(1:N))/1e6 M_seasonal_range(ind(1:N)) M_seasonal_error(ind(1:N)) M_max_month(ind(1:N)) M_min_month(ind(1:N))  ]

%%
col0 = hex2rgb('#CC2936'); 
col1 = hex2rgb('#4464AD'); 
col2 = hex2rgb('#F58F29'); 

figure('pos',[50 50 560*4.5/5.7  420*4.5/5.7])
ax=subplot(2,2,1) ;
%histogram((A_seasonal_range(seasonal))/1e6,0:.25:20,'facecolor',col0)
hst0=histogram(log10(A_seasonal_range(seasonal)/1e6),'facecolor',col0);
%title('Area','fontsize',7)
box off 
axis tight
xt   = 10.^(-3:3); 
set(gca,'fontsize',7,'xtick',log10(xt),'xticklabel',{'10^{-3}';'10^{-2}';'10^{-1}';'10^{0}';'10^{1}';'10^{2}';'10^{3}'}')
ntitle(' a ','fontsize',8,'fontweight','bold','location','nw')
text(-.15,40,{'Seasonal area range (km^2)'},'fontsize',7,'color',hst0(1).FaceColor,'horiz','center','vert','bot','fontangle','italic')
ylim([0 40])

ax(2)=subplot(2,2,2) ;
%histogram((M_seasonal_range(seasonal)),0:.25:15,'facecolor',col0)
hst0(2)=histogram(log10(M_seasonal_range(seasonal)),'facecolor',col0);
%title('Mass','fontsize',7)
box off 
axis tight
xt   = 10.^(-5:2); 
set(gca,'fontsize',7,'xtick',log10(xt),'xticklabel',{'10^{-5}';'10^{-4}';'10^{-3}';'10^{-2}';'10^{-1}';'10^{0}';'10^{1}';'10^{2}'}')
ntitle(' b ','fontsize',8,'fontweight','bold','location','nw')
text(-1.25,max(hst0(2).Values),{'Seasonal mass range (Gt)'},'fontsize',7,'color',hst0(2).FaceColor,'horiz','center','vert','bot','fontangle','italic')
ylim([0 40])

ax(3)=subplot(2,2,3);
hold on
hst(1) = histogram(A_max_month(seasonal),'facecolor',col1);
hst(2) = histogram(A_min_month(seasonal),'facecolor',col2);
axis([0.5 12.5 0 60])
set(gca,'fontsize',7,'xtick',1:12,'xticklabel',datestr(datenum(0,1:12,15),'mmm'))
ntitle(' c ','fontsize',8,'fontweight','bold','location','nw')


text(6,max(hst(1).Values),{'Maximum area'},'fontsize',7,'color',hst(1).FaceColor,'horiz','center','vert','bot','fontangle','italic')
text(10,max(hst(2).Values),{'Minimum area'},'fontsize',7,'color',hst(2).FaceColor,'horiz','center','vert','bot','fontangle','italic')

ax(4)=subplot(2,2,4);
hold on
hst(3)=histogram(M_max_month(seasonal),'facecolor',col1);
hst(4)=histogram(M_min_month(seasonal),'facecolor',col2);
axis([0.5 12.5 0 60])
set(gca,'fontsize',7,'xtick',1:12,'xticklabel',datestr(datenum(0,1:12,15),'mmm'))
ntitle(' d ','fontsize',8,'fontweight','bold','location','nw')

text(5,max(hst(3).Values),{'Maximum mass'},'fontsize',7,'color',hst(3).FaceColor,'horiz','center','vert','bot','fontangle','italic')
text(10,max(hst(4).Values),{'Minimum mass'},'fontsize',7,'color',hst(4).FaceColor,'horiz','center','vert','bot','fontangle','italic')

% 
ax(2).Position(1) = .5; 
ax(4).Position(1) = .5; 
ax(1).Position(2)=.5; 
ax(2).Position(2)=.5; 

% export_fig /Users/cgreene/Documents/GitHub/greenland-icemask/figures/ice_catchment_area_mass_seasonal_histograms.jpg -r600  -p0.01
% export_fig('/Users/cgreene/Documents/papers/greene2023greenland/figures/jpg/fig_ED03.jpg','-r900','-p0.01')
% exportgraphics(gcf,'/Users/cgreene/Documents/papers/greene2023greenland/figures/eps/fig_ED03.eps','ContentType','vector')

%% 

mo_max = nan(size(catchment),'single'); 

seas_range = mo_max; 
m_diff = mo_max; 

for k=1:261
    if has_data(k)
        m_diff(catchment==k & ice_sum==594) = M_diff_Gt(k); 
    end
    if seasonal(k)
        mo_max(catchment==k & ice_sum==594) = M_max_month(k); 
        seas_range(catchment==k & ice_sum==594) = M_seasonal_range(k); 
    end
end

alph = seas_range/4; 
alph(alph>1) = 1; 
figure
h=imagescn(x,y,mo_max);
h.AlphaData =  alph; 
bedmachine('gl','color',rgb('gray'),'linewidth',0.4,'greenland'); 
caxis([.5 12.5])
cmocean phase

readme = 'created by Chad Greene with ice_catchment_analysis_area_mass_timeseries.m, describing greenland_ice_masks_1972-2022_v1.nc'; 
%save('/Users/cgreene/Documents/GitHub/greenland-icemask/data/ice_catchment_analysis_area_mass_totals.mat','A_diff_km2','A_diff_km2_err','A_max_month','A_min_month','A_seasonal_range','names','lt','M_diff_Gt','M_diff_Gt_error','M_max_month','M_min_month','M_seasonal_range','has_data','readme')

%%

try 
    load('bedmachine_greenland_gl.mat')
catch
    [mask,xx,yy] = bedmachine_data('mask','greenland');
    [C,h] = contour(xx,yy,double(mask==1 | mask==2),[.5 .5],'k');
    
    [x_gl,y_gl] = C2xyz(C);
    x_gl = cell2nancat(x_gl); 
    y_gl = cell2nancat(y_gl); 
    %gl = polyshape(x1,y1);
    
    save('bedmachine_greenland_gl.mat','x_gl','y_gl')
end

%%*1.25

txtcol = 0.2*[1 1 1];

lowres = false; % use true for quick plots, false for full res

rockcol = 0.8*[1 1 1];

figure('pos',[150 500 720 292])
subsubplot(1,5,1)
if lowres
    h=imagescn(x,y,m_diff(1:4:end,1:4:end));
    maskoverlay(x(1:4:end),y(1:4:end),catchment(1:4:end,1:4:end)==0,'color',rockcol); 
else
    h=imagescn(x,y,m_diff);
    maskoverlay(x,y,catchment==0,'color',rockcol); 
end
hold on
plot(x_gl,y_gl,'color',rgb('gray'),'linewidth',0.25); 
axis image off  
cmocean -bal
caxis([-1 1]*150)
ntitle(' a ','location','nw','fontsize',8,'fontweight','bold')
cb = colorbar;
cb.Position = [.27 .2 .005 .2];
cb.FontSize = 6; 
cb(1).AxisLocation = 'in';
yl = ylabel(cb(1),{'Mass';'change';'1985-2022';'(Gt)'},'fontsize',6,'rotation',00);
yl(1).VerticalAlignment = 'middle'; 
yl(1).Position(1) = -6.5;
yl(1).Position(2) = -15;

txt(1) = text(-395042 +35e3 ,  -2312132+33e3,'Jakobshavn','fontangle','italic','fontsize',5,'horiz','right','color',txtcol,'vert','middle');
txt(2) = text(-466789+40e3 ,   -1020685,'Humboldt','fontangle','italic','fontsize',5,'horiz','right','color',txtcol);
txt(3) = text(566369  ,  -1078082+20e3,'Zachariae','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','bot');

plot([ -342703     -258940],[-2280583    -2276659],'-','color',txtcol,'linewidth',0.2)
plot([-423324     -395655],[-1037949    -1060214],'-','color',txtcol,'linewidth',0.2)
plot([505059      554138],[-1093016    -1056414+5e3],'-','color',txtcol,'linewidth',0.2)


subsubplot(1,5,2)
if lowres
    h=imagescn(x,y,seas_range(1:4:end,1:4:end));
    maskoverlay(x(1:4:end),y(1:4:end),catchment(1:4:end,1:4:end)==0,'color',rockcol); 
else
    h=imagescn(x,y,seas_range);
    maskoverlay(x,y,catchment==0,'color',rockcol); 
end
hold on
plot(x_gl,y_gl,'color',rgb('gray'),'linewidth',0.25); 
axis image off  
cmocean amp
%caxis([-1 1]*150)
ntitle(' b ','location','nw','fontsize',8,'fontweight','bold')
cb(2) = colorbar;
cb(2).Position = [.43 .2 .005 .2];
cb(2).FontSize = 6; 
cb(2).AxisLocation = 'in';
yl(2) = ylabel(cb(2),{'Seasonal';'mass';'range';'(Gt)'},'fontsize',6,'rotation',00);
yl(2).VerticalAlignment = 'middle'; 
yl(2).Position(1) = -6.5;
yl(2).Position(2) = 5.5;

%txt(4) = text(567963  ,  -2367963+30e3,'Kangerlussuaq','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','top');
%plot([ 508800      548986],[-2315051    -2337377],'-','color',txtcol,'linewidth',0.2)

txt(4) = text(567963  ,  -2153600-50e3,'Kangerlussuaq','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','bot');
plot([528579      559300],[-2202808    -2191490],'-','color',txtcol,'linewidth',0.2)

subsubplot(1,5,3)
if lowres
    h=imagescn(x,y,mo_max(1:4:end,1:4:end));
    maskoverlay(x(1:4:end),y(1:4:end),catchment(1:4:end,1:4:end)==0,'color',rockcol); 
    h.AlphaData =  alph(1:4:end,1:4:end); 
else
    h=imagescn(x,y,mo_max);
    maskoverlay(x,y,catchment==0,'color',rockcol); 
h.AlphaData =  alph; 
end
hold on
plot(x_gl,y_gl,'color',rgb('gray'),'linewidth',0.25); 
axis image off  
cmocean phase
caxis([0.5  12.5])
ntitle(' c ','location','nw','fontsize',8,'fontweight','bold')
axl = axis; 

addpath('/Users/cgreene/Documents/MATLAB/data_testing'); 
h_pb = itslive_phasebar('location','se','size',0.4,'color','w','centertext',{'Time of';'max extent'},'fontsize',5.5);
h_pb_ch = get(h_pb,'children'); 
set(h_pb_ch(1),'color','k')
h_pb.Position = [0.5330    0.14    0.07    0.2645]; 

good = M_diff_Gt<0 & M_seasonal_range>1e-5 & has_data' & ~lt & M_seasonal_range>M_seasonal_error & abs(A_diff_km2)>A_diff_km2_err;

%txt(5) = text(567963  ,  -2367963+10e3,'Helheimgletscher','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','top');



fp = plotboxpos(gca); 


foo=subplot(1,5,[4 5]);
scatter(M_seasonal_range(good),-M_diff_Gt(good),5,M_max_month(good),'filled')
%set(gca,'yaxislocation','right','xaxislocation','top')
caxis([0.5 12.5])
cmocean phase 
hold on

set(gca,'fontsize',7,'color','none','xscale','log','yscale','log')
%set(gca,'ydir','reverse')
xlabel('Range of seasonal mass variability (Gt)','fontsize',7)
ylabel('Mass loss from 1985 to 2022 (Gt)','fontsize',7)
axis tight
%axis([1e-5 14 .0007 520])
%axis([0.0014   14.0000    0.0171  297.3304])
set(gca,'xtick',10.^(-10:10),'ydir','reverse')

pv = polyfit(log10(M_seasonal_range(good)),log10(-M_diff_Gt(good)),1);
xtmp = linspace(min(M_seasonal_range(good)),max(M_seasonal_range(good)),1e3); 
pl = plot(xtmp,10.^(polyval(pv,log10(xtmp))),'--','color',.5*[1 1 1]); 
uistack(pl,'bottom')


axx = axes;
axx.InnerPosition = [0 0 1 1];
axis off
axx.Color = 'none';
text(fp(1)+fp(3)+.01,fp(2)+fp(4),' d ','horiz','left','vert','top','fontsize',8,'fontweight','bold')
uistack(axx,'bottom') 

foo.OuterPosition = [fp(1)+fp(3)+.01 fp(2) .28 fp(4)];

%export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/ice_catchment_mass_map_scatter.jpg','-r600','-p0.01')
% export_fig('/Users/cgreene/Documents/papers/greene2023greenland/figures/jpg/fig_03.jpg','-r900','-p0.01')
% exportgraphics(gcf,'/Users/cgreene/Documents/papers/greene2023greenland/figures/eps/fig_03.eps','ContentType','vector')
return

%% NO TEXT VERSION 


lowres = false; 

figure('pos',[150 500 720 292])
subsubplot(1,5,1)
if lowres
    h=imagescn(x,y,m_diff(1:4:end,1:4:end));
    maskoverlay(x(1:4:end),y(1:4:end),catchment(1:4:end,1:4:end)==0,'color',rockcol); 
else
    h=imagescn(x,y,m_diff);
    maskoverlay(x,y,catchment==0,'color',rockcol); 
end
hold on
plot(x_gl,y_gl,'color',rgb('gray'),'linewidth',0.25); 
axis image off  
cmocean -bal
caxis([-1 1]*150)
%ntitle(' a ','location','nw','fontsize',8,'fontweight','bold')
cb = colorbar;
cb.Position = [.27 .2 .005 .2];
cb.FontSize = 6; 
cb(1).AxisLocation = 'in';
yl = ylabel(cb(1),{'Mass';'change';'1985-2022';'(Gt)'},'fontsize',6,'rotation',00);
yl(1).VerticalAlignment = 'middle'; 
yl(1).Position(1) = -6.5;
yl(1).Position(2) = -15;


subsubplot(1,5,2)
if lowres
    h=imagescn(x,y,seas_range(1:4:end,1:4:end));
    maskoverlay(x(1:4:end),y(1:4:end),catchment(1:4:end,1:4:end)==0,'color',rockcol); 
else
    h=imagescn(x,y,seas_range);
    maskoverlay(x,y,catchment==0,'color',rockcol); 
end
hold on
plot(x_gl,y_gl,'color',rgb('gray'),'linewidth',0.25); 
axis image off  
cmocean amp
%caxis([-1 1]*150)
%ntitle(' b ','location','nw','fontsize',8,'fontweight','bold')
cb(2) = colorbar;
cb(2).Position = [.43 .2 .005 .2];
cb(2).FontSize = 6; 
cb(2).AxisLocation = 'in';
yl(2) = ylabel(cb(2),{'Seasonal';'mass';'range';'(Gt)'},'fontsize',6,'rotation',00);
yl(2).VerticalAlignment = 'middle'; 
yl(2).Position(1) = -6.5;
yl(2).Position(2) = 5.5;

%txt(4) = text(567963  ,  -2367963+30e3,'Kangerlussuaq','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','top');
%plot([ 508800      548986],[-2315051    -2337377],'-','color',txtcol,'linewidth',0.2)
% 
% txt(4) = text(567963  ,  -2153600-50e3,'Kangerlussuaq','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','bot');
% plot([528579      559300],[-2202808    -2191490],'-','color',txtcol,'linewidth',0.2)

subsubplot(1,5,3)
if lowres
    h=imagescn(x,y,mo_max(1:4:end,1:4:end));
    maskoverlay(x(1:4:end),y(1:4:end),catchment(1:4:end,1:4:end)==0,'color',rockcol); 
    h.AlphaData =  alph(1:4:end,1:4:end); 
else
    h=imagescn(x,y,mo_max);
    maskoverlay(x,y,catchment==0,'color',rockcol); 
h.AlphaData =  alph; 
end
hold on
plot(x_gl,y_gl,'color',rgb('gray'),'linewidth',0.25); 
axis image off  
cmocean phase
caxis([0.5  12.5])
%ntitle(' c ','location','nw','fontsize',8,'fontweight','bold')
axl = axis; 

addpath('/Users/cgreene/Documents/MATLAB/data_testing'); 
[h_pb,txtpb] = itslive_phasebar('location','se','size',0.4,'color','w','centertext',{'Time of';'max extent'},'fontsize',5.5);
h_pb_ch = get(h_pb,'children'); 
set(h_pb_ch(1),'color','k')
h_pb.Position = [0.5330    0.14    0.07    0.2645]; 

good = M_diff_Gt<0 & M_seasonal_range>1e-5 & has_data' & ~lt & M_seasonal_range>M_seasonal_error & abs(A_diff_km2)>A_diff_km2_err;

%txt(5) = text(567963  ,  -2367963+10e3,'Helheimgletscher','fontangle','italic','fontsize',5,'horiz','left','color',txtcol,'vert','top');

delete(txtpb(:))

fp = plotboxpos(gca); 


foo=subplot(1,5,[4 5]);
scatter(M_seasonal_range(good),-M_diff_Gt(good),5,M_max_month(good),'filled')
%set(gca,'yaxislocation','right','xaxislocation','top')
caxis([0.5 12.5])
cmocean phase 
hold on

set(gca,'fontsize',7,'color','none','xscale','log','yscale','log')
%set(gca,'ydir','reverse')
%xlabel('Range of seasonal mass variability (Gt)','fontsize',7)
%ylabel('Mass loss from 1985 to 2022 (Gt)','fontsize',7)
axis tight
%axis([1e-5 14 .0007 520])
%axis([0.0014   14.0000    0.0171  297.3304])
set(gca,'xtick',[],'ydir','reverse','ytick',[])

pv = polyfit(log10(M_seasonal_range(good)),log10(-M_diff_Gt(good)),1);
xtmp = linspace(min(M_seasonal_range(good)),max(M_seasonal_range(good)),1e3); 
pl = plot(xtmp,10.^(polyval(pv,log10(xtmp))),'--','color',.5*[1 1 1]); 
uistack(pl,'bottom')


axx = axes;
axx.InnerPosition = [0 0 1 1];
axis off
axx.Color = 'none';
%text(fp(1)+fp(3)+.01,fp(2)+fp(4),' d ','horiz','left','vert','top','fontsize',8,'fontweight','bold')
uistack(axx,'bottom') 

foo.OuterPosition = [fp(1)+fp(3)+.01 fp(2) .28 fp(4)];

cb(1).Ticks = []; 
cb(2).Ticks = []; 
cb(1).Label.String = ''; 
cb(2).Label.String = ''; 

% export_fig('/Users/cgreene/Documents/papers/greene2023greenland/figures/jpg/fig_03_NOTEXT.jpg','-r900','-p0.01')

%%

figure

seasonal = has_data' & A_seasonal_range>0*A_seasonal_error & M_seasonal_range>0*M_seasonal_error; 
sum(seasonal) 

scatter(M_seasonal_range(seasonal),-M_diff_Gt(seasonal),15,M_max_month(seasonal),'filled')
%set(gca,'yaxislocation','right','xaxislocation','top')
caxis([0.5 12.5])
cmocean phase 
hold on

set(gca,'fontsize',7,'color','none','xscale','log','yscale','log')

%%


% %yrs = 2014:2021;
% dA = nan(size(yrs)); 
% dM = dA; 
% 
% for k = 1:length(dA) 
%     dA(k) = A_tot(yr==yrs(k) & mo==m1) -  A_tot(yr==yrs(k) & mo==m2);
%     dM(k) = M_tot(yr==yrs(k) & mo==m1) -  M_tot(yr==yrs(k) & mo==m2);
% 
% end
% 
% figure
% plot(yrs,dM,'o-')

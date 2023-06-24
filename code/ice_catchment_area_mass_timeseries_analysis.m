
load('icemask_catchment_analysis_1972-2022_v1.mat') 

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
text(t1,tmp,[num2str(round(A_rate_km2_per_yr)),' km^{-2} yr^{-1}'],'horiz','left','vert','bot','color',col3,'fontsize',6)

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

%%

m1 = 5; % month of seasonal max 
m2 = 9; % month of seasonal min

[diff(As_median([m1 m2]))/1e6 hypot(As_std(m1),As_std(m2))/1e6 diff(Ms_median([m1 m2])) hypot(Ms_std(m1),Ms_std(m2))]


%%

%yrs = 2014:2021;
dA = nan(size(yrs)); 
dM = dA; 

for k = 1:length(dA) 
    dA(k) = A_tot(yr==yrs(k) & mo==m1) -  A_tot(yr==yrs(k) & mo==m2);
    dM(k) = M_tot(yr==yrs(k) & mo==m1) -  M_tot(yr==yrs(k) & mo==m2);

end

figure
plot(yrs,dM,'o-')

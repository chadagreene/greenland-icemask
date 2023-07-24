

load('ocean_data.mat','T_catchment_anomaly','T_catchment_anomaly_50');
load('sill_depth')
load('greenland_basin_runoff.mat'); 
load('ice_catchment_analysis_area_mass_totals.mat','M_diff_Gt','M_diff_Gt_error','M_max_month','M_seasonal_range');

load('testall.mat','Area','flux','sfz_slope','bed_changed','bed_slope','slope_mean','th_mean','tidewater','v_mean','width','x_deep','y_deep')

observation_data = ncread('greenland_ice_masks_1972-2022_v1.nc','observation_data');
has_data = (sum(observation_data>0)>5)';

S = shaperead('/Users/cgreene/Documents/data/basins/doi_10.7280_D1WT11__v1/Greenland_Basins_PS_v1.4.2.shp');

mt = strcmp({S.GL_TYPE},'TW')'; 
mt(261) = false; 

good = has_data & abs(M_diff_Gt)>M_diff_Gt_error & mt & M_diff_Gt<0 & width>0 &th_mean>0 & v_mean>0 & M_seasonal_range>2e-3  & isfinite(bed_slope) & isfinite(sill_depth); 

%%

good = good & isfinite(T_catchment_anomaly(:,10)); 

Number_of_catchments = sum(good) 

[CM,PM] = corrcoef([-log(-M_diff_Gt(good)) log(M_seasonal_range(good)) M_max_month(good) bed_slope(good) sfz_slope(good) bed_changed(good)   log(th_mean(good))     log(width(good))  log(runoff_basin(good)) log(v_mean(good))       log(flux(good))   sill_depth(good) T_catchment_anomaly(good,10)]); 
labM =             {'Mass change';           'Seasonal range';         'Seasonal timing'; 'Bed slope';   'Surface slope';'Bed elevation';   'Terminus thickness';   'Terminus width';   'Surface runoff';               'Terminus velocity';  'Terminus flux';  'Sill depth';    'T_{ocean}' }; 

CM_sq = round(CM.*abs(CM),2); 

% [CL,PL] = corrcoef([log(-M_diff_Gt(good)./Area(good)/917) log(M_seasonal_range(good)./Area(good)/917) M_max_month(good) bed_slope(good) sfz_slope(good) bed_changed(good) log(th_mean(good)) log(width(good)) log(Area(good)) log(v_mean(good)) log(flux(good)) sill_depth(good) log(runoff_basin_mmwe(good))]); 
% labL = {'dL'; 'L_{amp}'; 'L_{ph}'; 'bed slope'; 'surface slope';'bed elev'; 'thickness'; 'width';'area';'velocity';'flux'; 'sill depth'; 'runoff'}; 

[CL,PL] = corrcoef([-log(-M_diff_Gt(good)./Area(good)/917) log(M_seasonal_range(good)./Area(good)/917) M_max_month(good) bed_slope(good) sfz_slope(good) bed_changed(good)   log(th_mean(good))     log(width(good))  log(runoff_basin(good)) log(v_mean(good))       log(flux(good))   sill_depth(good) T_catchment_anomaly(good,10)]); 
labL =             {'Length change';                      'Seasonal range';                          'Seasonal timing'; 'Bed slope';   'Surface slope';'Bed elevation';   'Terminus thickness';   'Terminus width';   'Surface runoff';     'Terminus velocity';  'Terminus flux';    'Sill depth';    'T_{ocean}' }; 

CL_sq = round(CL.*abs(CL),2); 

%% 
warning off 

yl = [-160 -0.05];
hp= 0.02; % horiz pad between subplots
vp = 0.07; % vert pad between subplots 
mk = '.'; % maker
mks = 7; % markersize 
fs = 7; % fontsize 
fs2 = 7; 
axcol = .3*[1 1 1]; 

vpos = -350; 

%col = brewermap(12,'Dark2');
%col = hex2rgb({'#ee80fe', '#59a20c', '#762aac', '#6ab3e1', '#1c4c5e', '#48bea1', '#d10f55', '#096013', '#fe8f06', '#6c202e', '#da887c', '#413c09', '#b69cfd', '#1c4bb4', '#a7af87'});
col =  hex2rgb({'#782857','#68affc',  '#ab5031', '#387472', '#da9f63', '#4c319e', '#3bc185', '#de0ca3', '#454366', '#7fac0c', '#fc5468', '#056e12', '#cb84d4', '#6118df'});

figure('pos',[10 100 480*1.1 380*1.1])
subsubplot(3,4,1,'vpad',vp,'hpad',hp)
plot(M_seasonal_range(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(1,:))
axis tight 
ylim(yl)
set(gca,'xscale','log','yscale','log','xcolor',col(1,:),'ycolor',axcol,'fontsize',fs,'xtick',10.^(-10:10),'ytick',fliplr(-10.^(-10:10)))
fitline(M_seasonal_range(good),M_diff_Gt(good))
ntitle('Seasonal mass range','fontsize',fs,'fontweight','bold','vert','bot','color',col(1,:))
xlabel('Gt','fontsize',fs,'color',col(1,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,2)),'%.2f')];' '},'color',col(2-1,:),'fontsize',fs2,'location','sw')
ntitle({[' {\it p} = ',num2str(PM(1,2),1)]},'color',col(2-1,:),'fontsize',fs2,'location','sw')


subsubplot(3,4,2,'vpad',vp,'hpad',hp)
plot(M_max_month(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(2,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(2,:),'fontsize',fs,'xtick',1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'xticklabelrotation',0)
axis tight 
ylim(yl)
fitline(M_max_month(good),M_diff_Gt(good))
xlim([0.5 12.5])
%datetick('x','mmm','keeplimits')
ntitle('Seasonal timing','fontsize',fs,'fontweight','bold','vert','bot','color',col(2,:))
xlabel('Month of maximum mass','fontsize',fs,'color',col(2,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,3)),'%.2f')];' '},'color',col(3-1,:),'fontsize',fs2,'location','se')
ntitle({[' {\it p} = ',num2str(PM(1,3),1)]},'color',col(3-1,:),'fontsize',fs2,'location','se')

subsubplot(3,4,3,'vpad',vp,'hpad',hp)
plot(bed_slope(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(3,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(3,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(bed_slope(good),M_diff_Gt(good))
ntitle('Bed slope','fontsize',fs,'fontweight','bold','vert','bot','color',col(3,:))
xlabel('Degrees','fontsize',fs,'color',col(3,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,4)),'%.2f')];' '},'color',col(4-1,:),'fontsize',fs2,'location','se')
ntitle({[' {\it p} = ',num2str(PM(1,4),1)]},'color',col(4-1,:),'fontsize',fs2,'location','se')

subsubplot(3,4,4,'vpad',vp,'hpad',hp)
plot(sfz_slope(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(4,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(4,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(sfz_slope(good),M_diff_Gt(good))
ntitle('Surface slope','fontsize',fs,'fontweight','bold','vert','bot','color',col(4,:))
xlabel('Degrees','fontsize',fs,'color',col(4,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,5)),'%.2f')];' '},'color',col(5-1,:),'fontsize',fs2,'location','se')
ntitle({[' {\it p} = ',num2str(PM(1,5),1)]},'color',col(5-1,:),'fontsize',fs2,'location','se')

subsubplot(3,4,5,'vpad',vp,'hpad',hp)
plot(bed_changed(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(5,:))
set(gca,'xscale','lin','yscale','log','fontsize',fs,'xcolor',col(5,:),'ycolor',axcol,'ytick',fliplr(-10.^(-10:10)))
axis tight 
ylim(yl)
fitline(bed_changed(good),M_diff_Gt(good))
ntitle('Bed elevation','fontsize',fs,'fontweight','bold','vert','bot','color',col(5,:))
xlabel('Meters','fontsize',fs,'color',col(5,:))
ylabel('Total mass change 1985-2022 (Gt)','fontsize',fs); 
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,6)),'%.2f')];' '},'color',col(6-1,:),'fontsize',fs2,'location','se')
ntitle({[' {\it p} = ',num2str(PM(1,6),1)]},'color',col(6-1,:),'fontsize',fs2,'location','se')


subsubplot(3,4,6,'vpad',vp,'hpad',hp)
plot(th_mean(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(6,:))
set(gca,'xscale','log','yscale','log','xcolor',col(6,:),'fontsize',fs,'xtick',10.^(-10:10))
axis tight 
ylim(yl)
fitline(th_mean(good),M_diff_Gt(good))
ntitle('Terminus thickness','fontsize',fs,'fontweight','bold','vert','bot','color',col(6,:))
xlabel('Meters','fontsize',fs,'color',col(6,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,7)),'%.2f')];' '},'color',col(7-1,:),'fontsize',fs2,'location','sw')
ntitle({[' {\it p} = ',num2str(PM(1,7),1)]},'color',col(7-1,:),'fontsize',fs2,'location','sw')


subsubplot(3,4,7,'vpad',vp,'hpad',hp)
plot(width(good)/1000,M_diff_Gt(good),mk,'markersize',mks,'color',col(7,:))
set(gca,'xscale','log','yscale','log','xcolor',col(7,:),'fontsize',fs,'xtick',10.^(-10:10))
axis tight 
ylim(yl)
fitline(width(good)/1000,M_diff_Gt(good))
ntitle('Terminus width','fontsize',fs,'fontweight','bold','vert','bot','color',col(7,:))
xlabel('Kilometers','fontsize',fs,'color',col(7,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,8)),'%.2f')];' '},'color',col(8-1,:),'fontsize',fs2,'location','sw')
ntitle({[' {\it p} = ',num2str(PM(1,8),1)]},'color',col(8-1,:),'fontsize',fs2,'location','sw')


subsubplot(3,4,8,'vpad',vp,'hpad',hp)
plot(runoff_basin(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(8,:))
set(gca,'xscale','log','yscale','log','xcolor',col(8,:),'fontsize',fs,'xtick',10.^(-10:10))
axis tight 
ylim(yl)
fitline(runoff_basin(good),M_diff_Gt(good))
ntitle('Surface runoff','fontsize',fs,'fontweight','bold','vert','bot','color',col(8,:))
xlabel('Gt yr^{-1}','fontsize',fs,'color',col(8,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,9)),'%.2f')];' '},'color',col(9-1,:),'fontsize',fs2,'location','sw')
ntitle({[' {\it p} = ',num2str(PM(1,9),1)]},'color',col(9-1,:),'fontsize',fs2,'location','sw')


subsubplot(3,4,9,'vpad',vp,'hpad',hp)
plot(v_mean(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(9,:))
set(gca,'xscale','log','yscale','log','xcolor',col(9,:),'fontsize',fs,'ycolor',axcol,'xtick',10.^(-10:10),'ytick',fliplr(-10.^(-10:10)))
axis tight 
ylim(yl)
fitline(v_mean(good),M_diff_Gt(good))
ntitle('Terminus velocity','fontsize',fs,'fontweight','bold','vert','bot','color',col(9,:))
xlabel('m yr^{-1}','fontsize',fs,'color',col(9,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,10)),'%.2f')];' '},'color',col(10-1,:),'fontsize',fs2,'location','sw')
ntitle({[' {\it p} = ',num2str(PM(1,10),1)]},'color',col(10-1,:),'fontsize',fs2,'location','sw')

 

subsubplot(3,4,10,'vpad',vp,'hpad',hp)
plot(flux(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(10,:))
set(gca,'xscale','log','yscale','log','xcolor',col(10,:),'fontsize',fs,'xtick',10.^(-10:10))
axis tight 
ylim(yl)
fitline(flux(good),M_diff_Gt(good))
ntitle('Terminus flux','fontsize',fs,'fontweight','bold','vert','bot','color',col(10,:))
xlabel('Gt yr^{-1}','fontsize',fs,'color',col(10,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,11)),'%.2f')];' '},'color',col(11-1,:),'fontsize',fs2,'location','sw')
ntitle({[' {\it p} = ',num2str(PM(1,11),1)]},'color',col(11-1,:),'fontsize',fs2,'location','sw')


subsubplot(3,4,11,'vpad',vp,'hpad',hp)
plot(sill_depth(good),M_diff_Gt(good),mk,'markersize',mks,'color',col(11,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(11,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(sill_depth(good),M_diff_Gt(good))
ntitle('Sill depth','fontsize',fs,'fontweight','bold','vert','bot','color',col(11,:))
xlabel('Meters','fontsize',fs,'color',col(11,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,12)),'%.2f')];' '},'color',col(12-1,:),'fontsize',fs2,'location','se')
ntitle({[' {\it p} = ',num2str(PM(1,12),1)]},'color',col(12-1,:),'fontsize',fs2,'location','se')

subsubplot(3,4,12,'vpad',vp,'hpad',hp)
plot(T_catchment_anomaly(good,10),M_diff_Gt(good),mk,'markersize',mks,'color',col(12,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(12,:),'fontsize',fs,'xtick',-10:10)
axis tight 
ylim(yl)
fitline(T_catchment_anomaly(good,10),M_diff_Gt(good))
ntitle('Thermal forcing','fontsize',fs,'fontweight','bold','vert','bot','color',col(12,:))
xlabel('{\itT}_{ocean} ({\circ}C)','fontsize',fs,'color',col(12,:))
tmp = get(gca,'XLabel');
tmp.Position(2) = vpos;  
ntitle({[' {\it r}^2 = ',num2str(abs(CM_sq(1,13)),'%.2f')];' '},'color',col(13-1,:),'fontsize',fs2,'location','se')
ntitle({[' {\it p} = ',num2str(PM(1,13),1)]},'color',col(13-1,:),'fontsize',fs2,'location','se')

% export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_basin_mass_correlation_scatterplots.jpg','-pdf','-r600','-p0.01')

%%

figure('pos',[20 20 900 400]) 

subplot(1,2,1) 
set(gca,'outerpos',[0 0 .5 1]) 
%hm=heatmap(CM.*abs(CM));
hm=heatmap(CM);
caxis([-1 1])
crameri cork

%set(gca,'xtick',1:length(lab),'xticklabel',lab)
hm.XDisplayLabels = labM; 
hm.YDisplayLabels = labM; 
hm.CellLabelFormat = '%.2f';
hm.CellLabelColor = 1*[1 1 1];
hm.FontSize = 7; 

ax = gca; 
axp = struct(ax); 
axp.Axes.XAxisLocation = 'top'; 
colorbar off

subplot(1,2,2) 
set(gca,'outerpos',[0.46 0 .5 1]) 
%hm2=heatmap(CL.*abs(CL));
hm2=heatmap(CL);
caxis([-1 1])
crameri cork

%set(gca,'xtick',1:length(lab),'xticklabel',lab)
hm2.XDisplayLabels = labL; 
hm2.YDisplayLabels = labL; 
hm2.CellLabelFormat = '%.2f';
hm2.CellLabelColor = 1*[1 1 1];
hm2.FontSize = 7; 

ax2 = gca; 
axp2 = struct(ax2); 
axp2.Axes.XAxisLocation = 'top'; 
colorbar off
%colorbar


axx = axes('position',[0 0 1 1]);
axx.Color = 'none';
axis off
txt = text(0.25,.08,'Mass variability','horiz','center','fontsize',10,'fontweight','bold');
txt(2) = text(0.72,.08,'Effective length variability','horiz','center','fontsize',10,'fontweight','bold');

% export_fig /Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_basin_length_correlation_heatmap.jpg -pdf -r600 -p0.01

%%

function fitline(xx,yy)

tmpx = linspace(min(xx),max(xx),1000); 

if strcmpi(get(gca,'xscale'),'log')
    pv = polyfit(log10(xx),log10(yy),1);  
    tmpy  = 10.^(polyval(pv,log10(tmpx))); 
else
    pv = polyfit(xx,log10(yy),1);  
    tmpy  = 10.^(polyval(pv,tmpx)); 
end

hold on
plot(tmpx,tmpy,'-','color',.5*[1 1 1],'linewidth',0.4)
end

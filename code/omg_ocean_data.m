

%% 


fn = 'greenland_extruded_velocity_and_thickness_2023-04-06.nc';
x = ncread(fn,'x'); 
y = ncread(fn,'y'); 
catchment = permute(ncread(fn,'catchment'),[2 1]); 
rock = catchment==0; 
load('greenland_ice_masks_1972-2022_v1_ice_sum.mat')

[X,Y] = meshgrid(x,y); 
bed = bedmachine_interp('bed',X,Y,'greenland'); 

[mask_bm,x_bm,y_bm] = bedmachine_data('mask','greenland'); 

%%

perim = bwperim(rock | ice_sum>0) & ~rock;

influenced_ice = (ice_sum>0 & ice_sum<594) | perim; 

deepest_ice = nan(261,1);

for k = 1:261
    ind = influenced_ice & catchment==k; 
    if any(ind,'all')
        deepest_ice(k) = min(bed(ind)); 
    end

end

%%

fn = 'OMG_CTD_and_Float_Profiles.nc'; 

lat_omg = ncread(fn,'latitude'); 
lon_omg = ncread(fn,'longitude'); 
[x_omg,y_omg] = ll2psn(lat_omg,lon_omg); 
d_omg = -ncread(fn,'depths'); 
T_omg = ncread(fn,'potential_temperature')'; 
%S_omg = ncread(fn,'practical_salinity')'; 

T_omg(T_omg<-2.5) = NaN; % because there was a weird value around -4.5 that can't be true 
trim = d_omg<-995; % bc stuff gets weird below 1000 
T_omg(trim,:) = []; 
d_omg(trim) = []; 

basin_omg = interp2(x,y,catchment,x_omg,y_omg,'nearest',261);

dist2ice_km = interp2(x,y,bwdist(ice_sum>0),x_omg,y_omg)*0.120; 

T_max = max(T_omg,[],'omitnan'); 

T_bottom = nan(size(T_max)); 
d_bottom = nan(size(T_max)); 
T_deepest = nan(size(T_max)); 

T2 = T_omg; 

for k = 1:length(T_max)
   ind = find(isfinite(T_omg(:,k)),1,'last'); 
   d_bottom(k) = d_omg(ind); 
   T_bottom(k) = T_omg(ind,k); 

   
   basin_number = basin_omg(k);
   if basin_number>0
       T2(d_omg<deepest_ice(basin_number),k) = nan; 
       ind = find(isfinite(T2(:,k)),1,'last'); 
       if ~isempty(ind)
            T_deepest(k) = T2(ind,k); 
       end
   end
    
end

%%

cax = [-2 4];

T_mean_profile = movmean(mean(T_omg,2,'omitnan'),5); 

T_omg_anomaly = T_omg-T_mean_profile; 

T_mean_anomaly = mean(T_omg_anomaly,'omitnan'); 
T_mean_anomaly_50 = mean(T_omg_anomaly(d_omg<-50,:),'omitnan'); 

col = mat2rgb(T_mean_anomaly,cmocean('thermal'),[cax]); 

fs = 7; 
figure('pos',[100 500 546 298])

ax = subplot(1,2,1) ;
%plot(T_omg,d_omg,'linewidth',.25)
hold on
for k = 1:2828
    plot(T_omg(:,k),d_omg,'linewidth',.25,'color',col(k,:))
end
plot(T_mean_profile,d_omg,'k','linewidth',1)
box off
axis tight
ylim([-1000 0])
xlabel('Potential temperature (\circC)','fontsize',fs) 
ylabel('Depth (m)','fontsize',fs)
set(gca,'fontsize',fs) 

ax(2) = subplot(1,2,2); 
%bedmachine('gl','color',.5*[1 1 1],'greenland')

scatter(x_omg,y_omg,5,T_mean_anomaly,'filled')
hold on
axis tight off
axl = axis; 
mh(1) = maskoverlay(x_bm,y_bm,mask_bm==1,'color',.7*[1 1 1]);
hold on
mh(2) = maskoverlay(x_bm,y_bm,mask_bm==2,'color',.95*[1 1 1]);
uistack(mh(1),'bottom')
uistack(mh(2),'bottom'); 
axis(axl); 
ax(2).Position = [.45 .07 .33 .85];

cb = colorbar('southoutside'); 
cb.Position = [.54 .79 .15 .01];
%cb.Position = [.62 .4 .01 .2];
cb.FontSize = 6; 
cb.AxisLocation = 'in';
caxis([cax])
xlabel(cb,'Local thermal anomaly (K)')
cmocean thermal 

% cmocean bal 
% caxis([-6 6]) 
% set(cb,'xtick',-6:2:6)
% export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/omg_temperature_profiles_and_map.jpg','-r600','-p0.01') 
% export_fig('/Users/cgreene/Documents/papers/greene2023greenland/figures/jpg/fig_ED05.jpg','-r900','-p0.01')
% exportgraphics(gcf,'/Users/cgreene/Documents/papers/greene2023greenland/figures/eps/fig_ED05.eps','ContentType','vector')


%%


dist_threshold_km = 1:10; 
T_catchment_bottom = nan(261,length(dist_threshold_km)); 
T_catchment_deepest = T_catchment_bottom; 
T_catchment_anomaly = T_catchment_bottom; 
T_catchment_anomaly_50 = T_catchment_bottom; 

for kb = 1:261
    for kd = 1:length(dist_threshold_km) 
        ind = basin_omg==kb & dist2ice_km<=dist_threshold_km(kd);
        if any(ind)
            T_catchment_bottom(kb,kd) = max(T_bottom(ind));
            T_catchment_deepest(kb,kd) = max(T_deepest(ind));
            T_catchment_anomaly(kb,kd) = max(T_mean_anomaly(ind));
            T_catchment_anomaly_50(kb,kd) = max(T_mean_anomaly_50(ind));
        end
    end
end

readme = 'maximum bottom temperature observation corresponding to each glacier basin. Observations from NASA OMG. Created by omg_ocean_data.m'; 
%save('/Users/cgreene/Documents/GitHub/greenland-icemask/data/ocean_data.mat','dist_threshold_km','T_catchment_bottom','T_catchment_deepest','T_catchment_anomaly','T_catchment_anomaly_50','readme')

load('icemask_catchment_analysis_1972-2022_v1.mat') 

load('greenland_ice_masks_1972-2022_v1_ice_sum.mat')
load('greenland_icemask_length_analysis_v1.mat','outlets','xo','yo') 

obs = ncread('greenland_ice_masks_1972-2022_v1.nc','observation_data'); 
rock = permute(logical(ncread('greenland_ice_masks_1972-2022_v1.nc','rock')),[2 1]); 

fn = 'greenland_extruded_velocity_and_thickness_2023-04-06.nc';
vx = double(permute(ncread(fn,'vx'),[2 1])); 
vy = double(permute(ncread(fn,'vy'),[2 1])); 
th = permute(ncread(fn,'thickness'),[2 1]); 
th_err = permute(ncread(fn,'thickness_error'),[2 1]); 
catchment = permute(ncread(fn,'catchment'),[2 1]); 
v = hypot(double(vx),double(vy)); 

[X,Y] = meshgrid(x,y); 
bed = bedmachine_interp('bed',X,Y,'greenland'); 

%clear vx vy 

S = shaperead('/Users/cgreene/Documents/data/basins/doi_10.7280_D1WT11__v1/Greenland_Basins_PS_v1.4.2.shp');

%%
% 
% flux = nan(size(names)); 
% for k=1:261 
%     % x1 = S(k).X; 
%     % y1 = S(k).Y; 
%     % x1(end) = x1(1)+1; 
%     % y1(end) = y1(1)+1; 
%     
%     %[xi,yi] = mask2outline(x,y,imfill(never_ocean & catchment==k,'holes'));
%     
%     C = contourc(x,y,double(imfill(never_ocean & catchment==k,'holes')),[0.5 0.5]); 
% 
%     % Convert the contour matrix to x,y coordinates: 
%     [xb,yb] = C2xyz(C); 
% 
%     flux_tmp = nan(size(xb));
%     for kk = 1:length(xb)
% 
%     
%         [xi,yi] = psnpath([xb{kk} xb{kk}(1)+1],[yb{kk} yb{kk}(1)+1],10);
%         %[xi,yi] = psnpath(xi,yi,10);
%         
%         dxi = gradient(xi); 
%         dyi = gradient(yi); 
%         
%         hi  = interp2(x,y,th,xi,yi); 
%         vxi = interp2(x,y,vx,xi,yi); 
%         vyi = interp2(x,y,vy,xi,yi); 
%         flux_tmp(kk) = sum((vxi.*dyi -vyi.*dxi).*hi)*917*1e-12;
% 
%     end
%     flux(k) = sum(flux_tmp); 
%     k
% end

%% steady-state ice flux

rho_ice = 917;

% Determine the mass of ice in each grid cell: 
[Lat,~] = psn2ll(X,Y); 
A = (diff(x(1:2)) ./ psndistortion(Lat) ).^2; % grid cell area, (accounting for polar stereographic distortion).  
IceMass = rho_ice.*th.*A.*1e-12; % ice mass in each grid cell (Gt)
IceMass_plus = rho_ice.*(th + th_err).*A.*1e-12; % Thickness is primarily

clear Lat 

vx(rock) = 0; % should already be true, but just to be sure.
vy(rock) = 0; 

X1 = X; 
Y1 = Y; 
for k = 1:10
   dx = interp2(x,y,vx/10,X1,Y1); 
   dy = interp2(x,y,vy/10,X1,Y1); 
   
   X1 = X1+dx; 
   Y1 = Y1+dy; 
   k
end

never_ocean = ice_sum==594 & ismember(bedmachine_interp('mask',X,Y,'greenland'),[1 2]); 

% A grid of steady-state mass rate at the calving front: 
Mdot_cf = (double(never_ocean) - interp2(x,y,double(never_ocean),X1,Y1)) .* IceMass; 
Mdot_cf(isnan(Mdot_cf)) = 0; % bc we count it as zero anyway, and now we don't have to deal with nans. 

% % Mass rate with using v_plus: 
% Mdot_cf_plus_v = (double(never_ocean) - interp2(x,y,double(never_ocean),X1_plus,Y1_plus)) .* IceMass; 
% Mdot_cf_plus_v(isnan(Mdot_cf_plus_v)) = 0; 

% Mass rate using H_plus: 
Mdot_cf_plus_H = (double(never_ocean) - interp2(x,y,double(never_ocean),X1,Y1)) .* IceMass_plus; 
Mdot_cf_plus_H(isnan(Mdot_cf_plus_H)) = 0; 


flux = nan(size(names)); 
flux_plus = flux; 
for k = 1:261
    flux(k) = sum(Mdot_cf(never_ocean & catchment==k),'all');
    flux_plus(k) = sum(Mdot_cf_plus_H(never_ocean & catchment==k),'all');
    k
end

readme = 'created by greenland_basin_stats.m';
% save('greenland_basin_flux_steadystate.mat','flux','flux_plus','names','readme')
%%

% which basins have sufficient data for analysis? 
% requre at least 5 separate months throughout the record where terminus position was measured: 
has_data = sum(obs>0)>5; 

tidewater = strcmp({S.GL_TYPE},'TW')'; 
tidewater(261) = false; 

%% 

dx = diff(x(1:2));
width = nan(size(names)); 
th_mean = width; 
v_mean = width; 
slope_mean = width; 

for k = 1:261
    ind = outlets & catchment==k ; 
    
    width(k) = sum(ind,'all')*dx; 
    th_mean(k) = mean(th(ind)); 
    v_mean(k) = mean(v(ind)); 

end

%% 

outlet_deepest = nan(size(flux)); 
shallowest_fjord = nan(size(flux)); 

x_deep = nan(size(flux)); 
y_deep = nan(size(flux)); 

for k = 1:261
    tmp = v.*th; 
    tmp(~outlets | catchment~=k) = 0; 

    ind = find(tmp==max(tmp(:))); 
    if isscalar(ind)
        x_deep(k) = X(ind); 
        y_deep(k) = Y(ind); 
        outlet_deepest(k) = bed(ind); 
    end
    k
end

%%

bed_changed = nan(size(names)); 

changed = ice_sum>0 & ice_sum<594;

for k=1:261
    tmp = changed & catchment==k;
    if any(tmp,'all')
        bed_changed(k) = mean(bed(tmp)); 
    end
end

%%
% step = 0.3; % Number of streamline steps per grid cell.  
% reach = 500e3; % meters of reach in each direction from measured termini
% Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 
%        
% ds = stream2(x,y,vx,vy,x_deep,y_deep,[step Npts]);



%% 


step = 0.2; % Number of streamline steps per grid cell.  
reach = 7e3; % meters of reach in each direction from measured termini
Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 
       
ds = stream2(x,y,-vx,-vy,xo,yo,[step Npts]);

%%


x_5 = nan(size(ds))'; 
y_5 = nan(size(ds))'; 

for k = 1:length(ds) 
    
   xi = ds{k}(:,1); 
   yi = ds{k}(:,2); 
   if numel(xi)>1

        di = pathdistpsn(xi,yi);

      if di(end)>=5000
        x_5(k) = interp1(di,xi,5000); 
        y_5(k) = interp1(di,yi,5000); 
        
      end
   end
end

% Calculate straight-line distance to help identify crazies that could wrap back around too close to the end: 
d = hypot(x_5-xo,y_5-yo);


dbed = interp2(x,y,bed,x_5,y_5) - interp2(x,y,bed,xo,yo); 
good = d>4000 & d<6000 & isfinite(dbed); 
Dbed = gridbin(xo(good),yo(good),dbed(good),x,y,@mean); 

dsfz = interp2(x,y,bed+th,x_5,y_5) - interp2(x,y,bed+th,xo,yo); 
good = d>4000 & d<6000 & isfinite(dsfz); 
Dsfz = gridbin(xo(good),yo(good),dsfz(good),x,y,@mean); 

%% 

dbed_mean = nan(size(v_mean)); 
dsfz_mean = nan(size(v_mean)); 
Area = nan(size(v_mean)); 
for k = 1:261
    ind = outlets & catchment==k & isfinite(Dbed); 
    
    dbed_mean(k) = mean(Dbed(ind)); 
    dsfz_mean(k) = mean(Dsfz(ind)); 
    Area(k) = sum(th(ind))*120; 
end

%%

fn = '/Users/cgreene/Documents/data/SMB/gsfc_fdm_smb_v1_2_1_gris_Dec22.nc'; 
x_gsfc = ncread(fn,'x'); 
y_gsfc = ncread(fn,'y'); 
t_gsfc = ncread(fn,'time'); 
runoff = ncread(fn,'Ru');

s1 = squeeze(sum(sum(runoff,'omitnan'),'omitnan'))*grid_cell_area*rho_ice*1e-12;

t_range = t_gsfc>=2000 & t_gsfc<2021; 
mean_runoff_m_ice_per_year = mean(runoff,3,'omitnan'); % grid of mean rate of runoff.

grid_cell_area = 12500^2; % m^2
rho_ice = 917;            
kg_to_Gt = 1e-12; % multiplier for converting kg to Gt 
total_runoff_Gt_per_year = sum(mean_runoff_m_ice_per_year(:),'omitnan') * grid_cell_area * rho_ice * kg_to_Gt

t_range = t_gsfc>=2000 & t_gsfc<2021; 

runoff_annual = interp2(single(x_gsfc),single(y_gsfc),mean(runoff(:,:,t_range),3),X,Y) .*A * rho_ice * kg_to_Gt;

clear runoff t_gsfc x_gsfc y_gsfc fn

runoff_basin = nan(size(names)); 
for k = 1:261
    runoff_basin(k)  = sum(runoff_annual(catchment==k & ice_sum==594)); 
end

%%

M_ts_lp = movmean(M_ts,12); 
M_ts_hp = M_ts - M_ts_lp; 

ind_seas = t>=datenum(2013,9,15) & t<=datenum(2021,9,15); 

M_amp = nan(261,1); 
M_ph = M_amp; 
for k = 1:261
    ft = sinefit(t(ind_seas),M_ts_hp(ind_seas,k));
    M_amp(k) = ft(1); 
    M_ph(k) = ft(2); 
end

dM = ((M_ts(t==datenum(2021,9,15),:)-M_ts(t==datenum(1985,9,15),:)))'; 


%% 

% Glacier must have some width, some thickness, some motion, at least 5 months with terminus observations, and be nearly  
good  = width>0 &th_mean>0 & sum(obs>0)'>5 & v_mean>0 & isfinite(dbed_mean); 

tw = good & tidewater ; 
lt = good & ~tidewater; 


cla 
scatter(M_amp,-dM,30,M_ph,'filled')
caxis([1 365])
cmocean phase
set(gca,'xscale','log','yscale','log')
hold on
plot(M_amp(good),-dM(good),'ko')

%%

Z = bed + th;
Z(ice_sum<594) = bed(ice_sum<594);
Zg  = GRIDobj(x,y,Z);
Zgf  = fillsinks(Zg);

bed_filled = Zgf.Z; 

sill_depth  = nan(size(names)); 
for k = 1:261
    msk = catchment==k & ice_sum<594 & ice_sum>0; 
    if any(msk,'all')
        sill_depth(k) = min(bed_filled(msk)); 
    end
end

%%
basin_area = nan(size(names)); 
for k = 1:261
    basin_area(k) = sum(A(ice_sum==594 & catchment==k)); 
end

runoff_basin_mmwe = runoff_basin * 1e12 ./basin_area;

%%


bed_slope = atan2d(dbed_mean,5000); 
sfz_slope = atan2d(dsfz_mean,5000); 


good  = width>0 &th_mean>0 & sum(obs>0)'>5 & v_mean>0 & M_amp>2e-9  & isfinite(bed_slope); 

good = good & dM<0 ;

%C_log = corrcoef([log(-dM(good)./(th_mean(good).*width(good))) log(M_amp(good)./(th_mean(good).*width(good))) M_ph(good) dbed_mean(good) log(th_mean(good)) log(width(good)) log(v_mean(good)) log(th_mean(good).*width(good))  log(v_mean(good).*th_mean(good).*width(good))   ]); 


dL = dM./Area/917; 
L_amp = M_amp./Area/917; 

mass_display = false; 
if mass_display
    C = corrcoef([log(-dM(good)) log(M_amp(good)) M_ph(good) bed_slope(good) sfz_slope(good) bed_changed(good) log(th_mean(good)) log(width(good)) log(Area(good)) log(v_mean(good)) log(flux(good)) sill_depth(good) log(runoff_basin(good))]); 
    lab = {'dM'; 'M_{amp}'; 'M_{ph}'; 'bed slope'; 'surface slope';'bed elev'; 'thickness'; 'width';'area';'velocity';'flux'; 'sill depth'; 'runoff'}; 
else
    C = corrcoef([log(-dL(good)) log(L_amp(good)) M_ph(good) bed_slope(good) sfz_slope(good) bed_changed(good) log(th_mean(good)) log(width(good)) log(Area(good)) log(v_mean(good)) log(flux(good)) sill_depth(good) log(runoff_basin_mmwe(good))]); 
    lab = {'dL'; 'L_{amp}'; 'L_{ph}'; 'bed slope'; 'surface slope';'bed elev'; 'thickness'; 'width';'area';'velocity';'flux'; 'sill depth'; 'runoff'}; 

end

figure
clf
hm=heatmap(C);
caxis([-1 1])
crameri cork

%set(gca,'xtick',1:length(lab),'xticklabel',lab)
hm.XDisplayLabels = lab; 
hm.YDisplayLabels = lab; 
hm.CellLabelFormat = '%.2f';
hm.CellLabelColor = 1*[1 1 1];

ax = gca; 
axp = struct(ax); 
axp.Axes.XAxisLocation = 'top'; 

% export_fig greenland_basin_length_correlation_heatmap.jpg -r500 -p0.01


%%

yl = [-160 -0.0005];
hp= 0.02; % horiz pad between subplots
vp = 0.07; % vert pad between subplots 
mk = '.'; % maker
mks = 7; % markersize 
fs = 7; % fontsize 
axcol = .3*[1 1 1]; 

col = brewermap(12,'Dark2');

figure
subsubplot(3,4,1,'vpad',vp,'hpad',hp)
plot(M_amp(good),dM(good),mk,'markersize',mks,'color',col(1,:))
axis tight 
ylim(yl)
set(gca,'xscale','log','yscale','log','xcolor',col(1,:),'ycolor',axcol,'fontsize',fs)
fitline(M_amp(good),dM(good))
ntitle('M_{amp}','fontsize',fs,'fontweight','bold','vert','bot','color',col(1,:))
xlabel('Gt','fontsize',fs,'color',col(1,:))

subsubplot(3,4,2,'vpad',vp,'hpad',hp)
plot(M_ph(good),dM(good),mk,'markersize',mks,'color',col(2,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(2,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(M_ph(good),dM(good))
xlim([0 366])
datetick('x','mmm','keeplimits')
ntitle('M_{phase}','fontsize',fs,'fontweight','bold','vert','bot','color',col(2,:))
xlabel('day of max','fontsize',fs,'color',col(2,:))

subsubplot(3,4,3,'vpad',vp,'hpad',hp)
plot(bed_slope(good),dM(good),mk,'markersize',mks,'color',col(3,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(3,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(bed_slope(good),dM(good))
ntitle('bed slope','fontsize',fs,'fontweight','bold','vert','bot','color',col(3,:))
xlabel('degrees','fontsize',fs,'color',col(3,:))

subsubplot(3,4,4,'vpad',vp,'hpad',hp)
plot(sfz_slope(good),dM(good),mk,'markersize',mks,'color',col(4,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(4,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(sfz_slope(good),dM(good))
ntitle('surface slope','fontsize',fs,'fontweight','bold','vert','bot','color',col(4,:))
xlabel('degrees','fontsize',fs,'color',col(4,:))


subsubplot(3,4,5,'vpad',vp,'hpad',hp)
plot(bed_changed(good),dM(good),mk,'markersize',mks,'color',col(5,:))
set(gca,'xscale','lin','yscale','log','fontsize',fs,'xcolor',col(5,:),'ycolor',axcol)
axis tight 
ylim(yl)
fitline(bed_changed(good),dM(good))
ntitle('bed elevation','fontsize',fs,'fontweight','bold','vert','bot','color',col(5,:))
xlabel('meters','fontsize',fs,'color',col(5,:))

subsubplot(3,4,6,'vpad',vp,'hpad',hp)
plot(th_mean(good),dM(good),mk,'markersize',mks,'color',col(6,:))
set(gca,'xscale','log','yscale','log','xcolor',col(6,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(th_mean(good),dM(good))
ntitle('outlet thickness','fontsize',fs,'fontweight','bold','vert','bot','color',col(6,:))
xlabel('meters','fontsize',fs,'color',col(6,:))

subsubplot(3,4,7,'vpad',vp,'hpad',hp)
plot(width(good)/1000,dM(good),mk,'markersize',mks,'color',col(7,:))
set(gca,'xscale','log','yscale','log','xcolor',col(7,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(width(good)/1000,dM(good))
ntitle('outlet width','fontsize',fs,'fontweight','bold','vert','bot','color',col(7,:))
xlabel('kilometers','fontsize',fs,'color',col(7,:))

subsubplot(3,4,8,'vpad',vp,'hpad',hp)
plot(Area(good),dM(good),mk,'markersize',mks,'color',col(8,:))
set(gca,'xscale','log','yscale','log','xcolor',col(8,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(Area(good),dM(good))
ntitle('outlet area','fontsize',fs,'fontweight','bold','vert','bot','color',col(8,:))
xlabel('m^2','fontsize',fs,'color',col(8,:))

subsubplot(3,4,9,'vpad',vp,'hpad',hp)
plot(v_mean(good),dM(good),mk,'markersize',mks,'color',col(9,:))
set(gca,'xscale','log','yscale','log','xcolor',col(9,:),'fontsize',fs,'ycolor',axcol)
axis tight 
ylim(yl)
fitline(v_mean(good),dM(good))
ntitle('outlet velocity','fontsize',fs,'fontweight','bold','vert','bot','color',col(9,:))
xlabel('m/yr','fontsize',fs,'color',col(9,:))

subsubplot(3,4,10,'vpad',vp,'hpad',hp)
plot(flux(good),dM(good),mk,'markersize',mks,'color',col(10,:))
set(gca,'xscale','log','yscale','log','xcolor',col(10,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(flux(good),dM(good))
ntitle('glacier flux','fontsize',fs,'fontweight','bold','vert','bot','color',col(10,:))
xlabel('Gt/yr','fontsize',fs,'color',col(10,:))

subsubplot(3,4,11,'vpad',vp,'hpad',hp)
plot(sill_depth(good),dM(good),mk,'markersize',mks,'color',col(11,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(11,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(sill_depth(good),dM(good))
ntitle('sill depth','fontsize',fs,'fontweight','bold','vert','bot','color',col(11,:))
xlabel('meters','fontsize',fs,'color',col(11,:))

subsubplot(3,4,12,'vpad',vp,'hpad',hp)
plot(runoff_basin(good),dM(good),mk,'markersize',mks,'color',col(12,:))
set(gca,'xscale','log','yscale','log','xcolor',col(12,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(runoff_basin(good),dM(good))
ntitle('runoff','fontsize',fs,'fontweight','bold','vert','bot','color',col(12,:))
xlabel('Gt','fontsize',fs,'color',col(12,:))

%export_fig greenland_basin_mass_correlation_scatterplots.jpg -r500 -p0.01


%%

warning('off')

yl = [min(dL(good)) max(dL(good))];
hp= 0.02; % horiz pad between subplots
vp = 0.07; % vert pad between subplots 
mk = '.'; % maker
mkl = '.'; % land-terminating marker 
mks = 7; % markersize 
mksl = 7; 
fs = 7; % fontsize 
axcol = .3*[1 1 1]; 

col = brewermap(12,'Dark2');

figure
subsubplot(3,4,1,'vpad',vp,'hpad',hp)
plot(L_amp(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(1,:))
hold on
plot(L_amp(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(1,:))
axis tight 
ylim(yl)
set(gca,'xscale','log','yscale','log','xcolor',col(1,:),'ycolor',axcol,'fontsize',fs)
fitline(L_amp(good),dL(good))
ntitle('L_{amp}','fontsize',fs,'fontweight','bold','vert','bot','color',col(1,:))
xlabel('Gt','fontsize',fs,'color',col(1,:))

subsubplot(3,4,2,'vpad',vp,'hpad',hp)
plot(M_ph(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(2,:))
hold on
plot(M_ph(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(2,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(2,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(M_ph(good),dL(good))
xlim([0 366])
datetick('x','mmm','keeplimits')
ntitle('L_{phase}','fontsize',fs,'fontweight','bold','vert','bot','color',col(2,:))
xlabel('day of max','fontsize',fs,'color',col(2,:))

subsubplot(3,4,3,'vpad',vp,'hpad',hp)
plot(bed_slope(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(3,:))
hold on
plot(bed_slope(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(3,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(3,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(bed_slope(good),dL(good))
ntitle('bed slope','fontsize',fs,'fontweight','bold','vert','bot','color',col(3,:))
xlabel('degrees','fontsize',fs,'color',col(3,:))

subsubplot(3,4,4,'vpad',vp,'hpad',hp)
plot(sfz_slope(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(4,:))
hold on
plot(sfz_slope(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(4,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(4,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(sfz_slope(good),dL(good))
ntitle('surface slope','fontsize',fs,'fontweight','bold','vert','bot','color',col(4,:))
xlabel('degrees','fontsize',fs,'color',col(4,:))


subsubplot(3,4,5,'vpad',vp,'hpad',hp)
plot(bed_changed(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(5,:))
hold on
plot(bed_changed(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(5,:))
set(gca,'xscale','lin','yscale','log','fontsize',fs,'xcolor',col(5,:),'ycolor',axcol)
axis tight 
ylim(yl)
fitline(bed_changed(good),dL(good))
ntitle('bed elevation','fontsize',fs,'fontweight','bold','vert','bot','color',col(5,:))
xlabel('meters','fontsize',fs,'color',col(5,:))

subsubplot(3,4,6,'vpad',vp,'hpad',hp)
plot(th_mean(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(6,:))
hold on
plot(th_mean(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(6,:))
set(gca,'xscale','log','yscale','log','xcolor',col(6,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(th_mean(good),dL(good))
ntitle('outlet thickness','fontsize',fs,'fontweight','bold','vert','bot','color',col(6,:))
xlabel('meters','fontsize',fs,'color',col(6,:))

subsubplot(3,4,7,'vpad',vp,'hpad',hp)
plot(width(good & tidewater)/1000,dL(good & tidewater),mk,'markersize',mks,'color',col(7,:))
hold on
plot(width(good & ~tidewater)/1000,dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(7,:))
set(gca,'xscale','log','yscale','log','xcolor',col(7,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(width(good)/1000,dL(good))
ntitle('outlet width','fontsize',fs,'fontweight','bold','vert','bot','color',col(7,:))
xlabel('kilometers','fontsize',fs,'color',col(7,:))

subsubplot(3,4,8,'vpad',vp,'hpad',hp)
plot(Area(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(8,:))
hold on
plot(Area(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(8,:))
set(gca,'xscale','log','yscale','log','xcolor',col(8,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(Area(good),dL(good))
ntitle('outlet area','fontsize',fs,'fontweight','bold','vert','bot','color',col(8,:))
xlabel('m^2','fontsize',fs,'color',col(8,:))

subsubplot(3,4,9,'vpad',vp,'hpad',hp)
plot(v_mean(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(9,:))
hold on
plot(v_mean(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(9,:))
set(gca,'xscale','log','yscale','log','xcolor',col(9,:),'fontsize',fs,'ycolor',axcol)
axis tight 
ylim(yl)
fitline(v_mean(good),dL(good))
ntitle('outlet velocity','fontsize',fs,'fontweight','bold','vert','bot','color',col(9,:))
xlabel('m/yr','fontsize',fs,'color',col(9,:))

subsubplot(3,4,10,'vpad',vp,'hpad',hp)
plot(flux(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(10,:))
hold on
plot(flux(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(10,:))
set(gca,'xscale','log','yscale','log','xcolor',col(10,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(flux(good),dL(good))
ntitle('glacier flux','fontsize',fs,'fontweight','bold','vert','bot','color',col(10,:))
xlabel('Gt/yr','fontsize',fs,'color',col(10,:))
% 
subsubplot(3,4,11,'vpad',vp,'hpad',hp)
plot(sill_depth(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(11,:))
hold on
plot(sill_depth(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(11,:))
set(gca,'xscale','lin','yscale','log','xcolor',col(11,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(sill_depth(good),dL(good))
ntitle('sill depth','fontsize',fs,'fontweight','bold','vert','bot','color',col(11,:))
xlabel('meters','fontsize',fs,'color',col(11,:))
% 
subsubplot(3,4,12,'vpad',vp,'hpad',hp)
plot(runoff_basin_mmwe(good & tidewater),dL(good & tidewater),mk,'markersize',mks,'color',col(12,:))
hold on
plot(runoff_basin_mmwe(good & ~tidewater),dL(good & ~tidewater),mkl,'markersize',mksl,'color',col(12,:))
set(gca,'xscale','log','yscale','log','xcolor',col(12,:),'fontsize',fs)
axis tight 
ylim(yl)
fitline(runoff_basin_mmwe(good),dL(good))
ntitle('runoff','fontsize',fs,'fontweight','bold','vert','bot','color',col(12,:))
xlabel('mmwe','fontsize',fs,'color',col(12,:))

%export_fig greenland_basin_length_correlation_scatterplots.jpg -r500 -p0.01


%%

figure 
plot(bed_slope(good),log(-dM(good)),'x')
xlabel 'bed slope (degrees)'
hold on
plot(sfz_slope(good),log(-dM(good)),'o')

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
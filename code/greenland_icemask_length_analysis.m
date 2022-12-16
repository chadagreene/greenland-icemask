


%% Load data 
% Loading takes a few minutes 

fn = 'greenland_monthly_ice_masks_2022-11-18_clean2.nc';
ice = permute(logical(ncread(fn,'ice')),[2 1 3]); 
t = ncdateread(fn,'time'); 

filename = 'greenland_extruded_velocity_and_thickness_2022-11-17.nc'; 
x = double(ncread(fn,'x')); 
y = double(ncread(fn,'y')); 
vx = double(permute(ncread(filename,'vx'),[2 1])); 
vy = double(permute(ncread(filename,'vy'),[2 1])); 
th = double(permute(ncread(filename,'thickness'),[2 1])); 
rock = permute(ncread(filename,'v_source'),[2 1])==0; 

[X,Y] = meshgrid(x,y); 

%% Masks 
% This section takes 5 minutes 

always_ice = all(ice,3); 
never_ice = ~any(ice,3); 
sometimes_ice = ~always_ice & ~never_ice & ~rock; 

%% 

% Outlets are the ever-present edges of the ice sheet that are near active ice. 
outlets = bwperim(always_ice) & imdilate(sometimes_ice,strel('disk',3)); 

step = 0.1; % Number of streamline steps per grid cell.  
reach = 100e3; % meters of reach in each direction from measured termini
Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 
       
ds = stream2(x,y,vx,vy,X(outlets),Y(outlets),[step Npts]);

%%

dsc = ds; 
for k = 1:length(dsc) 
    dsc{k}(:,1) = k; 
    
   xi = ds{k}(:,1); 
   yi = ds{k}(:,2); 
   if numel(xi)>1
      dsc{k}(:,2) = pathdistpsn(xi,yi); 
      
   else 
       dsc{k}(:,2) = 0; 
   end
end

Ds = cell2mat(ds(:)); 
Dsc = cell2mat(dsc(:)); 

xi = Ds(:,1); 
yi = Ds(:,2); 
si = Dsc(:,1); % stream number
li = Dsc(:,2); 

%%

L = NaN(numel(ds),numel(t)); 

for kt=length(t):-1:1
    ind = interp2(x,y,ice(:,:,kt),Ds(:,1),Ds(:,2),'nearest'); 
    
    
    
    si_tmp = si(ind) ;
    li_tmp = li(ind); 

    parfor k = 1:numel(ds)
        %try
           % L(k,kt) = min(li(si==k & ~ind));
           % L(k,kt) = min(li_tmp(si_tmp==k)); 
           tmp = li_tmp(si_tmp==k);
           if ~isempty(tmp)
           L(k,kt) = tmp(end); 
           end
%         catch
           % L(k,kt) = NaN; 
        %end
    end
    disp([datestr(now),': Month ',num2str(kt),' of ',num2str(length(t))])
end

readme = 'Glacier length along flowlines starting at locations xo,yo. Created by greenland_icemask_length_analysis.m';
xo = X(outlets); 
yo = Y(outlets); 
%save('greenland_icemask_length_analysis.mat','xo','yo','t','L','readme')

%% Second try 

L = NaN(numel(ds),numel(t)); 

for kt=length(t):-1:1
    
    ind = interp2(x,y,ice(:,:,kt),Ds(:,1),Ds(:,2),'nearest'); 
    
    parfor k = 1:numel(ds)
        
        ind_tmp = ind(si==k);
        li_tmp = li(si==k); 
        ftmp = find(~[ind_tmp;false],1,'first'); 
        L(k,kt) = li_tmp(ftmp-1); 
        
    end
    disp([datestr(now),': Month ',num2str(kt),' of ',num2str(length(t))])
end

readme = 'Glacier length along flowlines starting at locations xo,yo. Created by greenland_icemask_length_analysis.m';
xo = X(outlets); 
yo = Y(outlets); 
%save('greenland_icemask_length_analysis3.mat','xo','yo','t','L','readme')

%%

L_lp = movmean(L,12,2,'omitnan');

dL = L_lp(:,t==datetime(2021,6,15)) - L_lp(:,t==datetime(1985,6,15));

ind = t>=datetime(2013,6,1) & t<=datetime(2021,6,15); 
A = nan(length(xo),1); 
ph = nan(length(xo),1); 
A_error = nan(length(xo),1); 


for k = 1:length(xo)

    indi = ind' & isfinite(L(k,:)) & isfinite(L_lp(k,:)); 
    
    ttmp = t(indi); 
    L_hp_tmp = L(k,indi) - L_lp(k,indi); 
        
%    ind2 = ~isoutlier(L_hp_tmp) & abs(L_hp_tmp)>0; 
    ind2 = abs(L_hp_tmp)>0; 
    
    if sum(ind2)>12
        ttmp = ttmp(ind2); 
        L_hp_tmp = L_hp_tmp(ind2); 

        ft = sinefit(ttmp,L_hp_tmp); 

        ind3 = ~isoutlier(L_hp_tmp - sineval(ft,ttmp')); 
        if sum(ind3)>12

            ttmp = ttmp(ind3); 
            L_hp_tmp = L_hp_tmp(ind3); 

            ft = sinefit(ttmp,L_hp_tmp); 

            A(k) = ft(1); 
            ph(k) = ft(2); 

            A_error(k) = std((L(k,indi) - L_lp(k,indi))' - sineval(ft,t(indi)));

        else
        A(k) = 0; 
        ph(k) = nan; 
        end
    else
        A(k) = 0; 
        ph(k) = nan; 
    end
 end

good = abs(dL)>0 & A>0 & A<10e3; 

figure
scatter(A(good)/1000,-dL(good)/1000,5,ph(good),'filled')
cmocean phase
caxis([0 365])
xlabel('seasonal length amplitude (km)')
ylabel('secular length change (km)') 

%%
tho = interp2(x,y,th,xo,yo); 

vo = interp2(x,y,hypot(vx,vy),xo,yo); 



step = 0.2; % Number of streamline steps per grid cell.  
reach = 1e6; % meters of reach in each direction from measured termini
Npts = round(reach/(diff(x(1:2))*step)); % Number of streamline points in each direction 

us = stream2(x,y,-vx,-vy,xo,yo,[step Npts]); 

% % Loop through each seed location to extrapolate terminus thickness along flowlines:  
% V = us; 
% S = us; 
% for k = 1:length(V) 
%    V{k}(:,1) = dL(k); 
%    S{k}(:,1) = A(k); 
%    S{k}(:,2) = ph(k); 
% end
% 
% % Concatenate the cells: 
% M = cell2mat(us(:)); % 
% V = cell2mat(V(:)); 
% 
% % Grid up the streamline data:  
% isf = isfinite(V(:,1)) & isfinite(M(:,1)); 
% Dlg = gridbin(M(isf,1),M(isf,2),V(isf,1),x,y); 
% 
% 
% holes = imfill(isfinite(Dlg),4,'holes') & ~isfinite(Dlg) & ~rock;
% Dlg = regionfill(Dlg,holes); 
% 
% 
% figure
% hdl = imagescn(x,y,Dlg/1000); 
% %hdl.AlphaData = hdl.AlphaData*.8; 
% axis image off
% cmocean -bal
% caxis([-1 1]*15)
% cb = colorbar; 
% ylabel(cb,'Glacier length change, 1985-2021 (km)')
% 
%%

col = mat2rgb(dL/1000,cmocean('-bal'),[-1 1]*15);

figure
hold on
for k = 3:5:length(dL)
    if good(k)
        plot(us{k}(:,1),us{k}(:,2),'color',col(k,:))
    end
    
end
box off
axis tight
daspect([1 1 1])

%%

col = mat2rgb(ph,cmocean('phase'),[0 365]);
%col = mat2rgb(A,cmocean('amp'),[0 2000]);

alph = A/1000; 
alph(alph>1) = 1; 

figure('pos',[10 10 800 1000])
hold on
for k = 1:1:length(dL)
    if true
        pl = plot(us{k}(:,1),us{k}(:,2),'color',col(k,:));
        pl.Color(4) = alph(k); 
    end
    
end
box off
axis tight off
daspect([1 1 1])
itslive_phasebar('location','se','color','w')

modismog('contrast','white')

% exportgraphics(gcf,'greenland_icemask_length_seasonality.jpg','resolution',600)

%%

close all
imagescn(x,y,hypot(vx,vy))
axis image 
hold on
plot(xo,yo,'k.')

%%

% Jakobshavn:
jako = inpolygon(xo,yo,[-180415.00    -180484.00    -180439.00    -180373.00],[-2277122.00   -2277149.00   -2277224.00   -2277182.00]);


% Petermann: 
petr = inpolygon(xo,yo,[-275754.00    -275747.00    -275666.00    -275666.00],[-949061.00    -949177.00    -949184.00    -949064.00]);

% Helheimgletcher: 
helh = inpolygon(xo,yo,[ 306393.00     306387.00     306442.00     306453.00],[-2577499.00   -2577551.00   -2577554.00   -2577503.00]);

% Zachariae: 
%zach = inpolygon(xo,yo,[ ],[]);

% Ryder: 
rydr = inpolygon(xo,yo,[-87670.00     -87729.00     -87638.00     -87612.00],[-885596.00    -885685.00    -885705.00    -885633.00]);

% 79 North? not sure if thats it
sevn = inpolygon(xo,yo,[464791.00     464732.00     464843.00     464883.00],[ -1021027.00   -1021178.00   -1021200.00   -1021093.00]);

djako = L(jako,t==datetime(1985,6,15))/1000 -0; 
dpetr = L(petr,t==datetime(1985,6,15))/1000 -9; 
dhelh = L(helh,t==datetime(1985,6,15))/1000 -30; 
%dzach = L(zach,t==datetime(1985,6,15))/1000 -30; 
drydr = L(rydr,t==datetime(1985,6,15))/1000 -40; 
dsevn = L(rydr,t==datetime(1985,6,15))/1000 -50; 


col_mean = hex2rgb('262730'); 
col_raw = hex2rgb('edcb96'); 
col_fit = hex2rgb('5762d5');

figure
plot(t,L(jako,:)/1000 - djako,'color',col_raw)
hold on
plot(t,L_lp(jako,:)/1000 - djako,'color',col_mean)
plot(t,L_lp(jako,:)/1000 + sineval([A(jako) ph(jako)],t')/1000 - djako,'color',col_fit)

plot(t,L(petr,:)/1000 - dpetr,'color',col_raw)
plot(t,L_lp(petr,:)/1000 - dpetr,'color',col_mean)
plot(t,L_lp(petr,:)/1000 + sineval([A(petr) ph(petr)],t')/1000 -dpetr,'color',col_fit)

plot(t,L(helh,:)/1000 - dhelh,'color',col_raw)
plot(t,L_lp(helh,:)/1000 - dhelh,'color',col_mean)
plot(t,L_lp(helh,:)/1000 + sineval([A(helh) ph(helh)],t')/1000 -dhelh,'color',col_fit)

plot(t,L(rydr,:)/1000 - drydr,'color',col_raw)
plot(t,L_lp(rydr,:)/1000 - drydr,'color',col_mean)
plot(t,L_lp(rydr,:)/1000 + sineval([A(rydr) ph(rydr)],t')/1000 -drydr,'color',col_fit)

plot(t,L(sevn,:)/1000 - dsevn,'color',col_raw)
plot(t,L_lp(sevn,:)/1000 - dsevn,'color',col_mean)
plot(t,L_lp(sevn,:)/1000 + sineval([A(sevn) ph(sevn)],t')/1000 -dsevn,'color',col_fit)

box off  
axis tight

% Patch indicating uncertain data: 
yl = ylim;
xl = xlim; 

hf = fill([xl(1) datetime(1985,6,15) datetime(1985,6,15) xl(1) xl(1)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'b');
hf.EdgeColor = 'none'; 
hf.FaceColor = 0.9*[1 1 1]; 
uistack(hf,'bottom')

hf = fill([datetime(2021,6,15) xl(2) xl(2) datetime(2021,6,15) datetime(2021,6,15)],[yl(1) yl(1) yl(2) yl(2) yl(1)],'b');
hf.EdgeColor = 'none'; 
hf.FaceColor = 0.9*[1 1 1]; 
uistack(hf,'bottom')

ylabel('Center flowline length (km)')

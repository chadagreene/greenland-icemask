% Run this after greenland_extrude_velocity_and_thickness.m.
%
% Run with R2021b because R2022b has issues with placing text objects. 
% 
% Chad A. Greene, NASA/JPL, November 2022. 


%%  Load data 

filename = 'greenland_extruded_velocity_and_thickness_2022-12-15.nc'; 

x = ncread(filename,'x'); 
y = ncread(filename,'y'); 

vx = double(permute(ncread(filename,'vx'),[2 1])); 
vy = double(permute(ncread(filename,'vy'),[2 1])); 
v_source = permute(ncread(filename,'v_source'),[2 1]); 
thickness = double(permute(ncread(filename,'thickness'),[2 1])); 
thickness_source = permute(ncread(filename,'thickness_source'),[2 1]); 
thickness_error = permute(ncread(filename,'thickness_error'),[2 1]); 
catchment = permute(ncread(filename,'catchment'),[2 1]); 

underlay = imfill(catchment<261,'holes'); 
catchment = double(catchment); 
catchment(catchment==0) = nan; 
catchment(catchment==261) = nan; 

%%

fst = 8; % fontsize, title 
fsl = 6; % fontsize, label

zl = [ -244567.60    -173694.31   -2334333.46   -2210298.49];

col = hex2rgb({'323031';'177e89';'5d2e8c';'db3a34';'ffc857'});
%vcol = hex2rgb({'323031';'177e89';'23C0D1';'5d2e8c';'db3a34';'ffc857'});

figure('pos',[20 50 669 556])
subsubplot(2,4,1) % subsubplot is a CDT function
imagescn(x,y,hypot(vx,vy))
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
set(gca,'colorscale','log')
caxis([1 10e3])
ax(1) = gca; 
cmocean thermal
ntitle('hypot(vx,vy)','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(1) = colorbar('south','color','w'); 
cb(1).Position(3:4) = cb(1).Position(3:4)*0.4; 
%cb(1).Position(1) = cb(1).Position(1)+.2; 
set(cb(1),'xtick',[1 10 100 1000 10000],'fontsize',fsl)
xlabel(cb(1),'ice speed (m yr^{-1})','fontsize',fsl)

plot(zl([1 2 2 1 1]),zl([3 3 4 4 3]),'linewidth',.3,'color','w')
%text(zl(2),zl(4),'e','horiz','left','vert','top','fontsize',fsl,'color','w')
ntitle(' a ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)

subsubplot(2,4,2)
imagescn(x,y,v_source)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
ax(2) = gca; 
ntitle('v\_source','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(2) = colorbar('south','color','k'); 
cb(2).Position(3:4) = cb(2).Position(3:4)*0.4; 
%cb(2).Position(1) = cb(2).Position(1)+.2; 
set(cb(2),'xtick',0:4,'fontsize',fsl,'xticklabel',{'rock','observed','interpolated','paleo flowlines','extrapolated'},'color','w')
colormap(gca,col)
caxis([-.5 4.5])
plot(zl([1 2 2 1 1]),zl([3 3 4 4 3]),'linewidth',.3,'color','w')
%text(zl(2),zl(4),'f','horiz','left','vert','top','fontsize',fsl,'color','w')
ntitle(' b ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)

subsubplot(2,4,5) 
imagescn(x,y,thickness)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
caxis([0 3300])
ax(3) = gca; 
cmocean ice
ntitle('thickness','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(3) = colorbar('south','color','w'); 
cb(3).Position(3:4) = cb(3).Position(3:4)*0.4; 
set(cb(3),'xtick',[0:1000:4000],'fontsize',fsl)
%cb(3).Position(1) = cb(3).Position(1)+.2; 
xlabel(cb(3),'ice thickness (m)','fontsize',fsl)
plot(zl([1 2 2 1 1]),zl([3 3 4 4 3]),'linewidth',.3,'color','w')
%text(zl(2),zl(4),'g','horiz','left','vert','top','fontsize',fsl,'color','w')
ntitle(' c ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)

subsubplot(2,4,6)
imagescn(x,y,thickness_source)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
ax(4) = gca; 
linkaxes(ax,'xy')
ntitle('thickness\_source','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(4) = colorbar('south','color','k'); 
cb(4).Position(3:4) = cb(4).Position(3:4)*0.4; 
set(cb(4),'xtick',[0:1000:4000],'fontsize',fsl)
%cb(4).Position(1) = cb(4).Position(1)+.2; 
set(cb(4),'xtick',0:3,'xticklabel',{'rock','BedMachine','AERODEM','inversion'},'color','w','fontsize',fsl)
colormap(gca,col([1 2 3 5],:))
caxis([-.5 3.5])
plot(zl([1 2 2 1 1]),zl([3 3 4 4 3]),'linewidth',.3,'color','w')
%text(zl(2),zl(4),'h','horiz','left','vert','top','fontsize',fsl,'color','w')
ntitle(' d ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)

[hs(1),ht(1)] = scalebarpsn('location','se','fontsize',fsl,'color','k');
hs(1).LineWidth = 1; 

subsubplot(2,4,3) 
imagescn(x,y,hypot(vx,vy))
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
set(gca,'colorscale','log')
caxis([1 15e3])
ax(5) = gca; 
cmocean thermal
ntitle('hypot(vx,vy)','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(5) = colorbar('south','color','w'); 
cb(5).Position(3:4) = cb(5).Position(3:4)*0.4; 
%cb(1).Position(1) = cb(1).Position(1)+.2; 
set(cb(5),'xtick',[1 10 100 1000 10000],'fontsize',fsl)
xlabel(cb(5),'ice speed (m yr^{-1})','fontsize',fsl)
ntitle(' e ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)
axis(zl) 

% Subset and plot velocity vectors: 
row = y>min(ylim) & y<max(ylim); 
c = x>min(xlim) & x<max(xlim); 
vx2 = vx(row,c); 
vy2 = vy(row,c); 
vs2 = v_source(row,c); 
x2 = x(c); 
y2 = y(row); 
hold on
q = quiversc(x2,y2,vx2,vy2,'density',100,'w');
q.AutoScaleFactor = 2; 
q.LineWidth = 0.35;


subsubplot(2,4,4)
imagescn(x,y,v_source)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
ax(6) = gca; 
ntitle('v\_source','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(6) = colorbar('south','color','k'); 
cb(6).Position(3:4) = cb(6).Position(3:4)*0.4; 
%cb(2).Position(1) = cb(2).Position(1)+.2; 
set(cb(6),'xtick',0:4,'fontsize',fsl,'xticklabel',{'rock','observed','interpolated','paleo flowlines','extrapolated'},'color','w')
colormap(gca,col)
caxis([-.5 4.5])
axis(zl)
ntitle(' f ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)

subsubplot(2,4,7) 
imagescn(x,y,thickness)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
caxis([0 1000])
ax(7) = gca; 
cmocean ice
ntitle('thickness','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(7) = colorbar('south','color','w'); 
cb(7).Position(3:4) = cb(7).Position(3:4)*0.4; 
set(cb(7),'xtick',[0:500:4000],'fontsize',fsl)
%cb(3).Position(1) = cb(3).Position(1)+.2; 
xlabel(cb(7),'ice thickness (m)','fontsize',fsl)
ntitle(' g ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)
axis(zl) 

subsubplot(2,4,8)
imagescn(x,y,thickness_source)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
ax(8) = gca; 
ntitle('thickness\_source','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(8) = colorbar('south','color','k'); 
cb(8).Position(3:4) = cb(8).Position(3:4)*0.4; 
set(cb(8),'xtick',[0:1000:4000],'fontsize',fsl)
%cb(4).Position(1) = cb(4).Position(1)+.2; 
set(cb(8),'xtick',0:3,'xticklabel',{'rock','BedMachine','AERODEM','inversion'},'color','w','fontsize',fsl)
colormap(gca,col([1 2 3 5],:))
caxis([-.5 3.5])
ntitle(' h ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)
axis(zl) 

[hs(2),ht(2)] = scalebarpsn('location','se','fontsize',fsl,'color','k');
hs(2).LineWidth = 1; 

set(gcf,'color','k')
% goal width 7.2 in = 183 mm

%export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_extruded_velocity_and_thickness_2022-12-15.jpg','-r600')

%%


figure('pos',[20 50 669 556])

subsubplot(2,4,5) 
imagescn(x,y,thickness_error)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
caxis([0 400])
ax = gca; 
cmocean amp
%ntitle('thickness error','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
cb(3) = colorbar('south','color','k'); 
cb(3).Position(3:4) = cb(3).Position(3:4)*0.4; 
set(cb(3),'xtick',[0:100:4000],'fontsize',fsl)
cb(3).Position(1)=.24;
%cb(3).Position(1) = cb(3).Position(1)+.2; 
xlabel(cb(3),'thickness error (m)','fontsize',fsl)
plot(zl([1 2 2 1 1]),zl([3 3 4 4 3]),'linewidth',.3,'color',rgb('blue'))
%text(zl(2),zl(4),'g','horiz','left','vert','top','fontsize',fsl,'color','w')

subsubplot(2,4,6) 
imagescn(x,y,thickness_error)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
caxis([0 400])
ax(2) = gca; 
cmocean amp
%ntitle('thickness','fontname','courier','color','w','fontweight','bold','fontsize',fst) 
%cb(7) = colorbar('south','color','w'); 
%cb(7).Position(3:4) = cb(7).Position(3:4)*0.4; 
%set(cb(7),'xtick',[0:500:4000],'fontsize',fsl)
%cb(3).Position(1) = cb(3).Position(1)+.2; 
%xlabel(cb(7),'ice thickness (m)','fontsize',fsl)
%ntitle(' g ','location','nw','fontsize',fst,'fontweight','bold','color','w','background','k','margin',0.001)
axis(zl) 

[hs(1),ht(1)] = scalebarpsn('location','se','fontsize',fsl,'color','k');
hs(1).LineWidth = 1; 

% export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_extruded_thickness_error_2022-12-15.jpg','-r600')

%%

rng('default')
N = 260; 
col = hsv2rgb([rand(N,1) 0.9*rand(N,1) 0.45+0.5*rand(N,1)]); 


figure('pos',[20 50 669 556])

subsubplot(2,4,5) 
hold on
maskoverlay(x,y,underlay,'color',.3*[1 1 1])
imagescn(x,y,catchment)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
caxis([.5 260.5])
colormap(gca,col)
plot(zl([1 2 2 1 1]),zl([3 3 4 4 3]),'linewidth',.3,'color',rgb('black'))

subsubplot(2,4,6) 
hold on
maskoverlay(x,y,underlay,'color',.3*[1 1 1])
imagescn(x,y,catchment)
axis image off
hold on
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
caxis([.5 260.5])
colormap(gca,col)
axis(zl) 

[hs(1),ht(1)] = scalebarpsn('location','se','fontsize',fsl,'color','k');
hs(1).LineWidth = 1; 

% export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_extruded_basins_2022-12-15.jpg','-r600')

% This script explores the uncertainty associated with inverting bed
% elevations for expected ice thickness. 
% 
% The script takes just a few seconds to run. 
% 
% Chad A. Greene, NASA/JPL, November 2022. 

%% Load data 

[th,x,y] = bedmachine_data('thickness','greenland'); 
bed = bedmachine_data('bed','greenland');  
errbed = bedmachine_data('errbed','greenland'); 
mask = bedmachine_data('mask','greenland'); 

%% Outlet locations 

% All pixels on the perimeter of the total ice sheet that aren't rock:  
outlets = bwperim(ismember(mask,[1 2])) & mask==2 & bed<0 & errbed<100; 

% Visual inspection: 
figure('pos',[63   596   599   268])
ax = subplot(1,2,1);

%axis([-218746     -184228    -2307014    -2245040])
axis([ -395936     -354905    -1537472    -1469087])
modismog('contrast','white')
hold on
maskoverlay(x,y,mask==1,'color',[0 0 0],'alpha',0.4)
maskoverlay(x,y,mask==2,'color',rgb('ice blue'),'alpha',0.6)
maskoverlay(x,y,outlets,'color',hex2rgb('#CB793A'))
q = itslive_quiver('region','GRE');
q.Color = hex2rgb('9A031E'); 
q.LineWidth = 0.3; 
axis off
[hsb(1),hsb(2)] = scalebarpsn('fontsize',5,'linewidth',1);
hsb(1).LineWidth = 1; 
%mapzoompsn('ne','frame','off')

%% Thickness inversion
% Ice at the balance point between grounded and floating has a thickness of
% -1.12*bed. 

th_equil = -1.12*bed; 

scale = median(th(outlets)./th_equil(outlets))

th_inv = th_equil*scale; 

inversion_error = std(th(outlets)-th_inv(outlets))


ax(2) = subplot(1,2,2);
histogram(th_equil(outlets)-th(outlets),-100:100,'edgecolor','none','facecolor',hex2rgb('#5F0F40'))
hold on
histogram(th_inv(outlets)-th(outlets),-100:100,'edgecolor','none','facecolor',hex2rgb('FCDC4D'))
box off
axis tight
xlabel('{\itH}_{inverted} - {\itH} (m)','fontsize',7)
ylabel('count','fontsize',7)
set(gca,'fontsize',7)
legend(['{\itH}_{eq.} (',num2str(median(th_equil(outlets)-th(outlets)),'%2.1f'),'±',num2str(std(th_equil(outlets)-th(outlets)),'%2.1f'),' m)'],...
  ['1.05*{\itH}_{eq.} (',num2str(median(th_inv(outlets)-th(outlets)),'%2.1f'),'±',num2str(std(th_inv(outlets)-th(outlets)),'%2.1f'),' m)'],...
  'fontsize',7,'location','northwest')
legend boxoff

ax(1).Position(1)=.25;

% exportgraphics(gcf,'/Users/cgreene/Documents/GitHub/greenland-coastlines/figures/thickness_inversion_uncertainty.jpg','resolution',600)

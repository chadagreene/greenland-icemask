% This script is just for quick tests to validate while writing 
% greenland_termini_to_icemasks.m. 
% 
% Chad A. Greene, July 2022. 

figure('color',.5*[1 1 1],'position',[42.00         55.00        560.00        892.00]); 
h = imagescn(x,y,land(:,:,1)); 
h.AlphaData = ~rock; 
axis image off
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
cmocean ice

txt = ntitle(datestr(double(t(1)),'yyyy'),'fontweight','bold','background','w'); 

[Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xlim,ylim,1000);
hold on
hh=image(xx,yy,Itmp(:,:,1:3));
uistack(hh,'bottom'); 
set(gca,'pos',[0 0 1 1])

scalebarpsn('location','se')
plup = plot(T.x(1),T.y(1),'b.'); % empty placeholder for now
pldn = plot(T.x(1),T.y(1),'y.'); % empty placeholder for now
%%

for k = 1:length(t)
   
   txt.String = datestr(double(t(k)),'yyyy');
   h.CData = land(:,:,k); 
   dt = T.t-t(k); 

   indup = dt>=0 & dt<30 & T.p>0;
   xup = T.x(indup) - interp2(x,y,vx,T.x(indup),T.y(indup)).*dt(indup)/365.25; 
   yup = T.y(indup) - interp2(x,y,vy,T.x(indup),T.y(indup)).*dt(indup)/365.25; 
   set(plup,'xdata',xup,'ydata',yup)


   inddn = dt<=0 & dt>-30 & T.p>0;
   xdn = T.x(inddn) - interp2(x,y,vx,T.x(inddn),T.y(inddn)).*dt(inddn)/365.25; 
   ydn = T.y(inddn) - interp2(x,y,vy,T.x(inddn),T.y(inddn)).*dt(inddn)/365.25; 
   set(pldn,'xdata',xdn,'ydata',ydn)
   
   drawnow
   pause(1/24)

end

%%
k=k+1; %1097 is insane.
% k=1052; 

txt.String = datestr(double(t(k)),'yyyy mm');
h.CData = land(:,:,k); 

dt = T.t-t(k); 

indup = dt>=0 & dt<30 & T.p>2;
xup = T.x(indup) - interp2(x,y,vx,T.x(indup),T.y(indup)).*dt(indup)/365.25; 
yup = T.y(indup) - interp2(x,y,vy,T.x(indup),T.y(indup)).*dt(indup)/365.25; 
set(plup,'xdata',xup,'ydata',yup)


inddn = dt<=0 & dt>-30 & T.p>2;
xdn = T.x(inddn) - interp2(x,y,vx,T.x(inddn),T.y(inddn)).*dt(inddn)/365.25; 
ydn = T.y(inddn) - interp2(x,y,vy,T.x(inddn),T.y(inddn)).*dt(inddn)/365.25; 
set(pldn,'xdata',xdn,'ydata',ydn)
drawnow

%%

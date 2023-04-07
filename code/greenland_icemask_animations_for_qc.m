


% Ice masks: 
fn = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_ice_masks_1972-2022_v0.nc';
x = double(ncread(fn,'x'));
y = double(ncread(fn,'y'));
t = ncdateread(fn,'time');
catchment = permute(ncread(fn,'catchment'),[2 1]); 

%t = datetime(1972,9,15):calmonths(1):datetime(2022,02,15);

tdn = datenum(t); 

D = load('terminus_data_densified_2023-01-09.mat');
color = parula(7); 

% catchi = zeros(size(D.x),'uint8'); 
% for k = 1:260
%     catchi(interp2(x,y,bwareafilt(catchment==k,1),D.x,D.y,'nearest')) = k; 
%     k
% end
load testcatch

%%
%for kd=278
for kd =256:260
   %clear A
   close all
   

   buf= 5e3; 

   try
   xl = [min(D.x(catchi==kd)) max(D.x(catchi==kd))] + buf*[-1 1];
   yl = [min(D.y(catchi==kd)) max(D.y(catchi==kd))] + buf*[-1 1];

   % Make the frame square: 
   dxl = diff(xl); 
   dyl = diff(yl); 
   if dxl>dyl
       yl = mean(yl) + dxl.*[-1 1]/2;
   else
       xl = mean(xl) + dyl.*[-1 1]/2; 
   end
   
   % Region of rows and columns of pixels to read: 
   ci=find((y>=yl(1))&(y<=yl(2)));
   ri=find((x>=xl(1))&(x<=xl(2)));
   
   x_tmp = x(ri); 
   y_tmp = y(ci); 
   
   ice_tmp = permute(logical(ncread(fn,'ice',[ri(1) ci(1) 1],[length(ri) length(ci) Inf])),[2 1 3]);
   rock_tmp = permute(logical(ncread(fn,'rock',[ri(1) ci(1)],[length(ri) length(ci)])),[2 1]);
   
   %
   
   figure('color',0*[1 1 1],'position',[28    87   931   782]); 
   h = imagescn(x_tmp,y_tmp,ice_tmp(:,:,1)); 
   h.AlphaData = ~rock_tmp; 
   axis image off
   bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
   cmocean ice
   
   txt = ntitle(num2str(year(t(1))),'fontweight','bold','backgroundcolor','w','fontsize',14); 
   
   [Itmp,xx,yy] = geoimread('/Users/cgreene/Documents/GreenlandBasemap50m.tif',xlim,ylim,1000);
   hold on
   hh=image(xx,yy,Itmp(:,:,1:3));
   uistack(hh,'bottom'); 
   set(gca,'pos',[0 0 1 1])
   
   scalebarpsn('location','se')
   pl1 = plot(mean(xl),mean(yl),'.','color',color(1,:)); % empty placeholder for now
   pl2 = plot(mean(xl),mean(yl),'.','color',color(2,:)); 
   pl3 = plot(mean(xl),mean(yl),'.','color',color(3,:)); 
   pl4 = plot(mean(xl),mean(yl),'.','color',color(4,:)); 
   pl5 = plot(mean(xl),mean(yl),'.','color',color(5,:)); 
   pl6 = plot(mean(xl),mean(yl),'.','color',color(6,:)); 
   pl7 = plot(mean(xl),mean(yl),'.','color',color(7,:)); 
   
   lg=label_greenland_glaciers('fontsize',10,'abbreviate');
   %
   xtxt = 1/14:1/7:1;
   for k=1:7
      lab(k) = text(xtxt(k),0,'','horiz','center','vert','bot','backgroundcolor',.5*[1 1 1],'color',color(k,:),'fontweight','bold','units','normalized','fontsize',14);
   end
   
   vf = VideoWriter(['/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/animations_quality_check/greenland_ice_masks_1972-2022_v1_',num2str(kd,'%03.f'),'.mp4'],'MPEG-4');
   vf.FrameRate = 12; 
   vf.Quality = 100; 
   open(vf)
   
   for k = 1:length(t)
      
      txt.String = datestr(t(k),'yyyy-mm-dd');
      h.CData = ice_tmp(:,:,k); 
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==1; 
      set(pl1,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(1).String = 'box'; 
      else
         lab(1).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==2; 
      set(pl2,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(2).String = 'manual'; 
      else
         lab(2).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==3; 
      set(pl3,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(3).String = 'slc-off'; 
      else
         lab(3).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==4; 
      set(pl4,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(4).String = 'auto id'; 
      else
         lab(4).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==5; 
      set(pl5,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(5).String = 'AutoTerm'; 
      else
         lab(5).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==6; 
      set(pl6,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(6).String = 'CALFIN'; 
      else
         lab(6).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==7; 
      set(pl7,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(7).String = 'MEaSUREs'; 
      else
         lab(7).String = ''; 
      end
   
   
      drawnow
      %pause(1/24)
   
      
      frame=getframe(gcf); 
      %frame = export_fig('-nocrop','-r200'); 
         writeVideo(vf,frame)
   
   end
   
   close(vf)
   end
end

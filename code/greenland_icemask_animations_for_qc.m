


% Ice masks: 
fn = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_monthly_ice_masks_2022-11-18.nc';
x = double(ncread(fn,'x'));
y = double(ncread(fn,'y'));
t = ncdateread(fn,'time');
tdn = datenum(t); 

D = load('terminus_data_densified_2022-11-14.mat');
color = parula(7); 

%%

for kd = 1:295
   clear A
   close all
   
   S_tmp = m_shaperead(['/Users/cgreene/Documents/data/coastlines/AutoTerm/GID',num2str(kd),'.shp']);
   t_tmp = S_tmp.dbfdata(:,1);
   
   for sk = 1:length(t_tmp)
   
      [A(sk,1).X,A(sk,1).Y] = ll2psn(S_tmp.ncst{sk}(:,2),S_tmp.ncst{sk}(:,1)); 
   
   end
   
   for k=1:length(A) 
      A(k).X = A(k).X'; 
      A(k).Y = A(k).Y'; 
   end
   
   
   
   buf = 5e3; 
   xl = [min([A.X]) max([A.X])] + buf*[-1 1];
   yl = [min([A.Y]) max([A.Y])] + buf*[-1 1];
   
   %
   
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
   
   
   vf = VideoWriter(['/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/animations_quality_check/greenland_icemask_synthesized_GID',num2str(kd),'_2022-11-18.mp4'],'MPEG-4');
   vf.FrameRate = 12; 
   vf.Quality = 100; 
   open(vf)
   
   
   
   for k = 1:length(t)
      
      txt.String = num2str(year(t(k)));
      h.CData = ice_tmp(:,:,k); 
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==1; 
      set(pl1,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(1).String = 'AutoTerm'; 
      else
         lab(1).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==2; 
      set(pl2,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(2).String = 'box'; 
      else
         lab(2).String = ''; 
      end
   
      ind = D.x>xl(1) & D.x<xl(2) & D.y>yl(1) & D.y<yl(2) & D.t>(tdn(k)-16) & D.t<(tdn(k)+16) & D.p==3; 
      set(pl3,'xdata',D.x(ind),'ydata',D.y(ind))
      if any(ind)
         lab(3).String = 'manual'; 
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
         lab(5).String = 'slc off'; 
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
         lab(7).String = 'Joughin'; 
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


load('icemask_catchment_analysis_1972-2022_v1.mat') 
names{244} = 'Ab Drachmann L Bistrup'; % just abbreviating for display purposes

%%
t = datetime(t,'convertfrom','datenum'); 
%ind = all(isfinite(A_ts),2) & t>=datetime(1985,8,15) & t<=datetime(2022,1,15); 
ind = true(size(t)); 

% seasonality may be best May 2015 thru May 2021

t = t(ind); 
A_ts = A_ts(ind,:); 
M_ts = M_ts(ind,:); 

%%
t=datenum(t); 


colt = hex2rgb('f9575a'); 
colc = hex2rgb('6f78f2'); 


lw2 = 1; % composite linewidth
ms=3; % markersize

figure('pos',[40 40 560 760])

for kk=0:6
   clf
   for k = 1:27

      ind_catchment = kk*27+k;
      if ind_catchment<=183
         subsubplot(7,4,k,'vpad',0.04,'hpad',0.04) 

         set(gca,'fontsize',5) 
            %boundedline(t,M_ts(:,ind_catchment),M_hi_ts(:,ind_catchment)-M_ts(:,ind_catchment))
            plot(t,M_ts(:,ind_catchment))
            hold on

            box off
            axis tight
            ntitle(names{ind_catchment},'fontsize',5,'color','k')

            ax = gca; 
            set(gca,'fontsize',5,'xcolor',.15*[1 1 1],'ycolor',.15*[1 1 1])

      end
      if kk==6
         k=21; % just to set the legend on the last plot properly
      end
   end

   % Make a legend for the last axes: 
   %lg = legend('response to thinning','response to calving','location','southwest');

   % Create a new axis just to get its position: 
   tmpax = subsubplot(7,4,k+1,'vpad',0.04,'hpad',0.04);

   % Move the legend to the empty axis position and delete the empty axes: 
   lg.Position=tmpax.Position;
   delete(tmpax)

   sgtitle({'Catchment mass (Gt)'},'fontsize',8) 

%   export_fig(['/Users/cgreene/Documents/GitHub/ice-shelf-geometry/figures/issm_timeseries/issm_GL_flux_from_thinning_and_calving_',num2str(kk+1),'.pdf'],'-r600','-painters','-p0.02')
   filename = '/Users/cgreene/Documents/GitHub/greenland-icemask/figures/test.pdf';
   if k==0
      assert(~exist(filename,'file'),[filename,' already exists'])
   end
   %export_fig(filename,'-r600','-painters','-nocrop','-append')

end

%%


M_sum = sum(M_ts-M_ts(1,:),2); 

A_sum = sum(A_ts-A_ts(1,:),2)/(1000^2); 

%%

figure('pos',[12.00        502.00        500.00        379.00]) 
hold on
set(gca,'fontsize',7)

Atmp = (M_ts-M_ts(1,:)); 

%[Atmp,ind] = sort(sum(Atmp,'omitnan')); 
[Atmp,ind] = sort(Atmp(end,:)); 
As = M_ts(:,ind) - M_ts(1,ind); 
isn = isnan(As(end,:)); 
As(end,isn) = As(end-1,isn); 

A_lo = fliplr(As(:,Atmp<0)); 
A_hi = As(:,Atmp>0); 
h_lo = area(t,A_lo,'edgecolor','none');
h_hi = area(t,A_hi,'edgecolor','none');

A_loc = [zeros(length(t),1) cumsum(A_lo,2)];
A_hic = [ zeros(length(t),1) cumsum(A_hi,2)];

name_2 = names(ind); 
name_lo = name_2(Atmp<0); 
name_lo = flipud(name_lo); 
name_hi = name_2(Atmp>0); 


cm = crameri('acton',22); 
cmind = 3+[1:5:19 2:5:19 3:5:19 4:5:19 5:5:19];
cm = repmat(cm(cmind,:),15,1); 
cm_lo = flipud(cm); 

for k = 1:length(h_lo)
   h_lo(k).FaceColor = cm(length(h_lo)-k+1,:); 
end
plot(t,A_loc(:,end),'color',cm(1,:),'linewidth',1); 

cm = crameri('oslo',25); 
cm(1:2,:) = []; 
cmind = 2+[1:5:19 2:5:19 3:5:19 4:5:19 5:5:19];
cm = repmat(cm(cmind,:),10,1); 
for k = 1:length(h_hi)
   h_hi(k).FaceColor = cm(k,:); 
end

%

Nlo = length(h_lo); 
A_loc = [zeros(length(t),1) cumsum(A_lo,2)]; 
yt = A_loc(end,2:end)-A_lo(end,:)/2; 
txt_lo=text(30+t(end)*ones(Nlo,1),yt(end-(Nlo-1):end),name_lo(end-(Nlo-1):end),...
    'horiz','left','vert','middle','fontsize',6,'fontangle','italic','color',.15*[1 1 1]);

box off
axis tight


shadowcol = [0 0 0];
linecol = hex2rgb('f0d079'); 
plot(t,M_sum,'color',shadowcol,'linewidth',1.1)
plot(t(end),M_sum(end),'.','linewidth',1.1,'color',shadowcol,'markersize',8)
plot(t(1),0,'.','linewidth',1.1,'color',shadowcol,'markersize',8)
plot(t,M_sum,'color',linecol,'linewidth',.8)
plot(t(end),M_sum(end),'.','linewidth',.8,'color',linecol,'markersize',7)
plot(t(1),0,'.','linewidth',.8,'color',linecol,'markersize',7)

txtsiz = linspace(0.1,6,Nlo); 
txtsiz = [linspace(0.1,2,150),linspace(2.1,6,Nlo-150)];
for k = 1:Nlo
    txt_lo(k).FontSize = txtsiz(k); 
    %txt_lo(k).Color = cm_lo(k,:); 
end
ylim([-1060 5])
text(t(end)+50,M_sum(end),'Greenland ','color',linecol,'fontname','Helvetica','fontsize',7,'fontweight','bold','horiz','left','vert','mid')

%
ax = axis; 
axcol = 0.6*[1 1 1]; 
xl = xlim; 
yl = ylim; 
hl = plot(xl,repmat((-2000:100:200)',1,2),'-','color',axcol); 
uistack(hl,'bottom')
txt = text(xl(1)*ones(length(-2000:100:200),1),(-2000:100:200)',num2str((-2000:100:200)'),'horiz','left','vert','bot','fontsize',7,'color',axcol,'clipping','on');
uistack(txt(:),'bottom')
%txt(end).String = '+20\times1000 km^2'; 
%txt(end-1).String = '+10'; 
%txt(22).String = '+100 Gt'; 
txt(21).String = ''; 
txt(11).String = '-1000 Gt';
text(xl(1),yl(end),'Cumulative mass change due to glacier retreat since 1985','fontsize',7,'horiz','left','vert','bot','fontweight','bold','color',.1*[1 1 1])
set(gca,'ycolor','none','xcolor',.15*[1 1 1])
axis(ax)
xtick = datenum(1986:2:2022,1,1); 
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,'yyyy'))

gc = gca; 
gc.Position(1) = 0.02;

% export_fig('/Users/cgreene/Documents/GitHub/greenland-icemask/figures/greenland_cumulative_masschange_1985-2021.jpg','-r900','-p0.01')

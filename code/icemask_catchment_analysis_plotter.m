
load('icemask_catchment_analysis_clean.mat') 
names{244} = 'Ab Drachmann L Bistrup'; % just abbreviating for display purposes

%% Seasonal 

M_ts_lp = movmean(M_ts,12); 
M_ts_hp = M_ts - M_ts_lp; 


%%
% 
ind = all(isfinite(A_ts),2) & t>=datetime(1985,6,15) & t<=datetime(2021,6,15); 

% ind = 154:586; % for highest confidence data June 1985 through June 2021. 

t = t(ind); 
A_ts = A_ts(ind,:); 
M_ts = M_ts(ind,:); 
t=datenum(t); 

M_ts_lp = M_ts_lp(ind,:); 
M_ts_hp = M_ts_hp(ind,:); 

%%


M_amp = nan(1,261); 
M_ph = M_amp; 
for k = 1:261
    ft = sinefit(t,M_ts_hp(:,k));
    M_amp(k) = ft(1); 
    M_ph(k) = ft(2); 
end

dM = (M_ts(end,:)-M_ts(1,:)); 

figure
scatter(M_amp,dM,30,M_ph,'filled')
caxis([0 365])
cmocean phase 

%%

% figure
% plot(t,M_ts-M_ts(1,:))
% hold on
% plot(t,sum(M_ts-M_ts(1,:),2),'k','linewidth',2)


M_sum = sum(M_ts-M_ts(1,:),2); 

A_sum = sum(A_ts-A_ts(1,:),2)/(1000^2); 
% 
% figure
% plot((A_ts(end,:)-A_ts(1,:))/1000^2,M_ts(end,:)-M_ts(1,:),'.')
% box off
% axis tight
% text((A_ts(end,:)-A_ts(1,:))/1000^2,M_ts(end,:)-M_ts(1,:),names,'fontsize',7,'vert','bot','horiz','center'); 

%%
if false 
dM = (M_ts(end,:)-M_ts(1,:)); 

f0 = find(dM==0); 

extruded_filename = 'greenland_extruded_velocity_and_thickness_2022-11-29.nc'; 
icemask_filename = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_monthly_ice_masks_2022-11-18.nc';

x = double(ncread(extruded_filename,'x')); 
y = double(ncread(extruded_filename,'y')); 

rock = permute(ncread(extruded_filename,'v_source'),[2 1])==0; 
catchment = double(permute(ncread(extruded_filename,'catchment'),[2 1])); 

ice = permute(logical(ncread(icemask_filename,'ice',[1 1 500],[Inf Inf 1])),[2 1]);

fn_term = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/terminus_data_densified_2022-11-14.mat';
T = load(fn_term); 

figure
imagescn(x,y,ice & ismember(catchment,f0))
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
hold on
plot(T.x(1:3:end),T.y(1:3:end),'r.')
end
%%
if false
dM_map = nan(size(ice)); 
for k = 1:261
    dM_map(catchment==k & ice) = dM(k); 
end

figure
imagescn(x,y,dM_map)
bedmachine('gl','color',rgb('gray'),'linewidth',.3,'greenland')
hold on
colorbar
caxis([-1 1]*max(abs(dM)))
cmocean -bal
end
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

textborder(t(end)+50,M_sum(end),'Greenland ',linecol,shadowcol,'fontname','Helvetica','fontsize',7,'fontweight','bold','horiz','left','vert','mid')

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
txt(22).String = '+100 Gt'; 
%txt(11).String = ''; 
text(xl(1),yl(end),'Cumulative mass change due to glacier retreat since 1985','fontsize',7,'horiz','left','vert','bot','fontweight','bold','color',.1*[1 1 1])
set(gca,'ycolor','none','xcolor',.15*[1 1 1])
axis(ax)
xtick = datenum(1985:2:2021,6,15); 
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,'yyyy'))

gc = gca; 
gc.Position(1) = 0.02;





%% SUBFUNCTION 

function h = textborder(x, y, string, text_color, border_color, varargin)
%TEXTBORDER Display text with border.
%   TEXTBORDER(X, Y, STRING)
%   Creates text on the current figure with a one-pixel border around it.
%   The default colors are white text on a black border, which provides
%   high contrast in most situations.
%   
%   TEXTBORDER(X, Y, STRING, TEXT_COLOR, BORDER_COLOR)
%   Optional TEXT_COLOR and BORDER_COLOR specify the colors to be used.
%   
%   Optional properties for the native TEXT function (such as 'FontSize')
%   can be supplied after all the other parameters.
%   Since usually the units of the parent axes are not pixels, resizing it
%   may subtly change the border of the text out of position. Either set
%   the right size for the figure before calling TEXTBORDER, or always
%   redraw the figure after resizing it.
%   
%   Author: João F. Henriques, April 2010
% Heavily edited by Chad Greene for FOP tools, November 2015, including the
% following changes: 
%   * Now prints 8 background texts rather than the previous 4. 
%   * Now returns object handles 
%   * Can now handle multiple text inputs
%   * Now prints in data units rather than data units, converting to pixels, moving, 
%       then converting back to data units. This is b/c changing units was VERY slow
%       on R2014b for some reason.  Took about 35 seconds before, which is absurd.  



if nargin < 5, border_color = 'k'; end  %default: black border
if nargin < 4, text_color = 'w'; end  %default: white text

pos = getpixelposition(gca); 
xl = get(gca,'xlim'); 
xperpx = diff(xl)/pos(3); 
yl = get(gca,'ylim'); 
yperpx = diff(yl)/pos(4); 

% border around the text, composed of 8 text objects  
%offsets = [xperpx yperpx].*[[0 -1; -1 0; 0 1; 1 0] ; 0.71*[ 1 1; -1 1; -1 -1; 1 -1]];
offsets = [xperpx yperpx].*[[0 -1; -1 0; 0 1; 1 0] ; .9*[ 1 1; -1 1; -1 -1; 1 -1]];

offsets = offsets/3;
% Initialize counters for background and main text:  
cbg = 1; 
    
for n = 1:length(x)     
	for k = 1:8
		h.bg(cbg) = text(x(n)+offsets(k,1), y(n)+offsets(k,2), string, 'Color',border_color, varargin{:});
		
        % Increment background counter: 
        cbg = cbg+1; 
   end

end

% the actual text inside the border
h.t = text(x, y, string, 'Color',text_color, varargin{:});

if nargout==0
    clear h;
end

end
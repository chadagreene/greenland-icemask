% A quick estimate of how much ice loss has occurred above or below the 
% surface of hydrostatic equilibrium. 

%%

fn = 'greenland_extruded_velocity_and_thickness_2023-04-06.nc';
x = ncread(fn,'x');
y = ncread(fn,'y');
catchment = permute(ncread(fn,'catchment'),[2 1]);

th = permute(ncread(fn,'thickness'),[2 1]);
th_error = permute(ncread(fn,'thickness_error'),[2 1]);
load('greenland_ice_masks_1972-2022_v1_ice_sum.mat')

[X,Y] = meshgrid(x,y);
bed = bedmachine_interp('bed',X,Y,'greenland');
sfz = bed+th;

%%
freeboard = base2freeboard(bed); % base2freeboard is a function in AMT
freeboard(bed>=0) = bed(bed>=0); % Accounts for bedrock above sea level.

taf = th - (freeboard-bed);

taf_high = th + th_error - (freeboard-bed);

%

changed = ice_sum>0 & ice_sum<594 & sfz>0;

% Neglecting polar stereographic distorion and just taking the max minus
% min extents of the whole ice sheet but this is just for a minor bit of
% context in the paper so it should be adequate for this case. Neglecting
% ps distortion will cause way less error than the error estimates for this
% particular 

MassLoss_Total = sum(th(changed))*120*120*917*1e-12

MassLoss_Above_Flotation = sum(taf(changed & taf>0))*120*120*917*1e-12

%% 

Percent_measurable = MassLoss_Above_Flotation*100/1034


%%

changed & taf>0

%% The following is not a great way to estimate uncertainty 
% the estimates below are too high so maybe don't use them. 

massloss_error = nan(261,1) ; 
for k=1:261
    
    massloss_error(k) = sum(taf_high(changed & taf_high>0 & catchment==k))*120*120*917*1e-12 - sum(taf(changed & taf>0 & catchment==k))*120*120*917*1e-12; 

end

MassLoss_Above_Flotation_error = rssq(massloss_error)

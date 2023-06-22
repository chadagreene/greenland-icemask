
load('icemask_catchment_analysis_1972-2022_v1.mat') 

A_err = A_err'; 
A_ts = A_ts'; 
M_err_pick_ts = M_err_pick_ts'; 
M_err_th_ts = M_err_th_ts'; 
M_ts = M_ts'; 

M_ts_err = hypot(M_err_th_ts,M_err_pick_ts); 


A_interann = movmean(A_ts,12,2); 
A_seasonal = A_ts - A_interann; 

M_interann = movmean(M_ts,12,2); 
M_seasonal = M_ts - M_interann; 

%t = datetime(t,'convertfrom','datenum'); 

M_tot = sum(M_ts); 
M_tot_err = rssq(M_ts_err); 

A_tot = sum(A_ts); 
A_tot_err = rssq(A_err); 
figure

boundedline(datenum(t),A_tot,A_tot_err,'nan','gap')

figure
boundedline(datenum(t),M_tot,M_tot_err,'nan','gap')
hold on
plot(datenum(t),sum(M_interann))

%% 

tmp = sum(M_seasonal,2); 

[yr,mo,dy] = datevec(t); 

yrs = 2014:2020; 
col = parula(length(yrs)); 

Ms = nan(12,length(yrs)); 

figure
hold on
for k = 1:12

    ind = yr>=yrs(1) & yr<=yrs(end) & mo==k; 
    plot(k,tmp(ind),'.','color',.4*[1 1 1])

    Ms(k,:) = tmp(ind); 
end

[hb1,hb2] = boundedline(1:12,median(Ms,2),std(Ms,[],2)); 
uistack(hb1,'bottom'); 
uistack(hb2,'bottom'); 

%%
tmp = median(Ms,2);

tmp(9)-tmp(5)


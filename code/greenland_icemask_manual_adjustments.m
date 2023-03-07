

% Ice masks: 
fn = '/Users/cgreene/Documents/data/coastlines/greenland-coastlines-greene/greenland_monthly_ice_masks_2023-02-22.nc';
x = double(ncread(fn,'x'));
y = double(ncread(fn,'y'));
t = ncdateread(fn,'time');
t = datetime(1972,9,15):calmonths(1):datetime(2022,02,15);

tdn = datenum(t); 

D = load('terminus_data_densified_2023-01-09.mat');
color = parula(7); 
%%
kd = 150
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
   
   
   
   buf = 9e3; 
   xl = [min([A.X]) max([A.X])] + buf*[-1 1];
   yl = [min([A.Y]) max([A.Y])] + buf*[-1 1];
   
   %
   
   % Region of rows and columns of pixels to read: 
   ci=find((y>=yl(1))&(y<=yl(2)));
   ri=find((x>=xl(1))&(x<=xl(2)));
   
   x_tmp = x(ri); 
   y_tmp = y(ci); 
   
   ice = permute(logical(ncread(fn,'ice',[ri(1) ci(1) 1],[length(ri) length(ci) Inf])),[2 1 3]);
   rock = permute(logical(ncread(fn,'rock',[ri(1) ci(1)],[length(ri) length(ci)])),[2 1]);



sometimes_ice = any(ice,3) & ~all(ice,3); 
ice_sum = sum(ice(:,:,154:586),3); 
[X,Y] = meshgrid(x_tmp,y_tmp); 

flicker = sum(abs(double(ice(:,:,1:end-1))-double(ice(:,:,2:end))),3);

figure
subplot(1,2,1) 
imagescn(x_tmp,y_tmp,ice_sum)
axis image 
hold on
ax = gca; 
subplot(1,2,2) 
imagescn(x_tmp,y_tmp,flicker)
axis image 
ax(2)  = gca; 
linkaxes(ax,'xy'); 

%% 2 

trim = sometimes_ice & X>-263217 & X<-257850 & Y>-1984775 & Y<-1979589 & ice_sum<100; 

for k = 154:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 9

trim = sometimes_ice & X>-318146 & X<-314490 & Y>-1826849 & Y<-1823349 & ice_sum<86 ;
for k = 1:239
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 10

trim = sometimes_ice & ice_sum<433 & inpolygon(X,Y,[-322183     -321592     -322678     -323954     -323935],[-1807295    -1805943    -1805638    -1805314    -1806762]); 

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 27 

trim = ~ice(:,:,248) & sometimes_ice & ice_sum<433 & inpolygon(X,Y,[-367902     -364253     -358703     -359615     -362884],[ -1578989    -1582790    -1581270    -1575796    -1575872]); 

for k = 248:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

trim = ~ice(:,:,331) & ~ice(:,:,333) & sometimes_ice & ice_sum<433 & inpolygon(X,Y,[-367902     -364253     -358703     -359615     -362884],[ -1578989    -1582790    -1581270    -1575796    -1575872]); 

for k = 331:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 36

trim = sometimes_ice & ice_sum<300 & inpolygon(X,Y,[-378634     -373851     -373375     -374803     -377777],[-1496243    -1496410    -1494173    -1492745    -1492817]); 

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 37
trim = sometimes_ice & ice_sum<200 & inpolygon(X,Y,[-385219     -382166     -380802     -380750     -382803     -385271],[-1482959    -1482907    -1481335    -1479503    -1479218    -1479672]); 

for k = 1:384
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

trim = sometimes_ice & ice_sum<290 & ~ice(:,:,385) & inpolygon(X,Y,[-385219     -382166     -380802     -380750     -382803     -385271],[-1482959    -1482907    -1481335    -1479503    -1479218    -1479672]); 

for k = 385:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 39 tongue 

trim = sometimes_ice & ice_sum<100 & inpolygon(X,Y,[ -388348     -388950     -389837     -390361     -391191     -390224     -387745],[  -1476823    -1477028    -1477619    -1478233    -1479223    -1479928    -1476971]); 

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 40 

trim = sometimes_ice & ice_sum<230 & ~ice(:,:,122) & inpolygon(X,Y,[ -390876     -394700     -396164     -393480     -388273],[ -1462088    -1465206    -1466888    -1469301    -1468732]); 
trim = trim | (sometimes_ice & ice_sum<300& inpolygon(X,Y,[-395475     -395486     -394811     -394764],[ -1466838    -1467496    -1467543    -1466826 ])); 

for k = 122:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 48

trim = sometimes_ice & ice_sum<260 & ~ice(:,:,380) & inpolygon(X,Y,[-439894     -439527     -436162     -434052     -436162 ],[-1423335    -1432053    -1431471    -1426118    -1422448]); 

for k = 380:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 53 total failure 

disp 'need to address 53'

%% 60 

filler = sometimes_ice & inpolygon(X,Y,[ -524269     -524245     -522356     -521734     -520921     -523719],[-1403707    -1404544    -1408060    -1408395    -1408012    -1403803]); 

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(filler) = true; 
    ice(:,:,k) = tmp; 
end

%% 63

filler = sometimes_ice & inpolygon(X,Y,[-553586     -553608     -552135     -551920     -552957],[ -1383762    -1384026    -1384059    -1383629    -1383657]);

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(filler) = true; 
    ice(:,:,k) = tmp; 
end

%% 64

filler = sometimes_ice & inpolygon(X,Y,[-568553	-568115	-567360	-567754	-567841	-567885	-568039	-568280	-568838	-568849],[-1386506	-1386418	-1386298	-1385772	-1385640	-1385137	-1384874	-1384622	-1384709	-1385936]);

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(filler) = true; 
    ice(:,:,k) = tmp; 
end

%% 65
% some flicker remains in later years. 

trim = sometimes_ice & ice_sum<590 & inpolygon(X,Y,[ -570454     -571455     -571788     -570909     -566725     -568847     -569484     -569848],[-1391790    -1392882    -1394034    -1395004    -1390698    -1389364    -1390031    -1391365]);

for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 68

trim = sometimes_ice & ice_sum<590 & inpolygon(X,Y,[-569886     -570042     -565293     -565395     -565476     -565788     -565881],[ -1331144    -1331579    -1332287    -1331434    -1331068    -1330565    -1329173]); 
for k = 1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 70

trim = sometimes_ice & ice_sum<590 & inpolygon(X,Y,[-591095     -588758     -589531     -589550     -589801     -589743     -591153],[-1287004    -1287120    -1286522    -1285672    -1285150    -1284610    -1285092]); 
trim = trim | (ice_sum<10 & inpolygon(X,Y,[ -591230     -588527     -588082     -591481],[ -1287333    -1287275    -1283663    -1283683])); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 84

trim = sometimes_ice & ~ice(:,:,71) & inpolygon(X,Y,[ -599056     -594851     -594956     -599612],[  -1139474    -1138570    -1136034    -1136763]); 

for k=1:70
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 85 Humboldt 

trim = sometimes_ice & ice_sum<310 & inpolygon(X,Y,[ -398812     -398754     -411675     -411501     -408125],[-1086615    -1082424    -1087837    -1091970    -1096917]); 
trim =  trim | (sometimes_ice & ice_sum<400 & inpolygon(X,Y,[ -403992     -403294     -398870     -398870],[-1083646    -1089059    -1086964    -1082482]));

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 97 

trim = sometimes_ice & ice_sum<260 & inpolygon(X,Y,[188301	188852	189402	189677	189585	189524	189096	188301	191756	193895	196402	196127	192520	190319],[-869117	-869759	-870768	-871655	-872480	-873152	-873703	-874406	-875445	-875903	-871960	-869056	-869331	-868475]); 

trim = trim | inpolygon(X,Y,[187829      189301      195055      195055      191523],[ -869089     -870079     -870882     -867697     -867322]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 108 

%trim = sometimes_ice &  inpolygon(X,Y,[],[]); 
trim = sometimes_ice &  inpolygon(X,Y,[  604352      604527      605295      605723      605681],[-1437733    -1438629    -1438649    -1437984    -1437616]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 121 

trim = sometimes_ice &  inpolygon(X,Y,[586407      586035      585553      585954      587119      586718],[  -2056048    -2056572    -2057342    -2057738    -2056166    -2055898]); 
for k=1:157
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

trim = ice_sum>0 &  inpolygon(X,Y,[588946      588642      588981      589527      589335],[-2060549    -2061322    -2061691    -2060683    -2060541]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 123 

trim = ice_sum>0 &  inpolygon(X,Y,[605339	605363	606319	606845	608160	611315	614853	610551	608949	608566	608423	608088	607562	606702	606224],[-1971757	-1972498	-1974124	-1975128	-1976323	-1979264	-1974267	-1969463	-1970036	-1970299	-1970801	-1971255	-1971662	-1972020	-1971829]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 127 

trim = ice_sum>0 &  inpolygon(X,Y,[729451	729458	729780	730460	731118	731170	730729	729825],[-2026156	-2027397	-2028331	-2028787	-2028697	-2028271	-2026567	-2026111]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end
%% 128 

trim = ice_sum>0 &  inpolygon(X,Y,[781668	781968	782340	782640	782981	783291	783581	784005	784118	784222	784263	784294	784511	784449	783084],[-2025499	-2026398	-2026305	-2026419	-2026688	-2026915	-2026967	-2026626	-2026347	-2026078	-2025757	-2025592	-2025396	-2024320	-2024351]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 138
trim = ice_sum>0 &  inpolygon(X,Y,[708446	708391	710397	710807	710160	709732	709604	709713	709558	709294	709157	709057	708875	708711	708610],[-2231889	-2232071	-2233530	-2232992	-2230613	-2230485	-2230604	-2230786	-2231014	-2231142	-2231342	-2231388	-2231670	-2231716	-2231853]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 139
trim = ice_sum>0 &  inpolygon(X,Y,[704979      704491      704159      703758      703496      703741      705136      705363],[  -2233012    -2233099    -2233134    -2233309    -2233518    -2234216    -2233727    -2233221]); 

for k=188:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

trim = ice_sum>0 &  inpolygon(X,Y,[704811	704690	704359	704004	703665	703584	703382	703213	703067	703415	703826	705110	705029],[-2232576	-2232665	-2232714	-2232819	-2232875	-2233053	-2233247	-2233327	-2233465	-2234232	-2234167	-2233255	-2232722]); 

for k=298:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 141
trim = ice_sum>0& sometimes_ice &  inpolygon(X,Y,[ 668695      668729      669911      669658      670633      670094      669670      669498      669372],[-2261073    -2263665    -2263722    -2262243    -2261578    -2261268    -2261142    -2261027    -2260844]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%% 150
trim = ice_sum>0&  inpolygon(X,Y,[504721	504477	504404	504319	504197	504209	504099	503965	504258	505662	508054	507896	507896	507847	507725	507542	507359	507139	507053	507224	507383	507554	507627	507700	507591	507444	507273	507066	506907	504892],[-2285430	-2285504	-2285662	-2285833	-2285992	-2286285	-2286456	-2286578	-2287091	-2287750	-2287933	-2287726	-2287579	-2287481	-2287323	-2287262	-2287176	-2287079	-2286822	-2286651	-2286358	-2286126	-2285968	-2285797	-2285650	-2285577	-2285565	-2285382	-2285125	-2285235]); 

for k=1:length(t)
    tmp = ice(:,:,k); 
    tmp(trim) = false; 
    ice(:,:,k) = tmp; 
end

%%


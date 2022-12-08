
% This is a quick script for visually inspecting AutoTerm data and setting
% error thresholds 
% 
% Chad Greene, November 2022. 

d = dir('/Users/cgreene/Documents/data/coastlines/AutoTerm/*.shp');

%%
figure
clear A ta sensor err
   %S_tmp = m_shaperead(['/Users/cgreene/Documents/data/coastlines/AutoTerm/',d(kd).name]);
   S_tmp = m_shaperead('/Users/cgreene/Documents/data/coastlines/AutoTerm/GID86.shp'); 
   t_tmp = S_tmp.dbfdata(:,1);
   s_tmp = S_tmp.dbfdata(:,3);
   e_tmp = S_tmp.dbfdata(:,4);


   for sk = 1:length(t_tmp)
   
      [A(sk,1).X,A(sk,1).Y] = ll2psn(S_tmp.ncst{sk}(:,2),S_tmp.ncst{sk}(:,1)); 
   
      ta(sk,1) = datenum(t_tmp{sk});
      %sensor{sk,1} = s_tmp{sk}; 
      err(sk,1) = str2double(e_tmp{sk}); 
   end

   for k=1:length(A) 
   A(k).X = A(k).X'; 
   A(k).Y = A(k).Y'; 
end
cla

ue = unique(err)
col = parula(length(ue)); 
hold on
cla
for k=length(ue):-1:1
   plot([A(err==ue(k)).X],[A(err==ue(k)).Y],'.','color',col(k,:),'markersize',1)
end
%title(num2str(kd))
axis tight
daspect([1 1 1])
axis([axis + 50e3*[-1 1 -1 1]])
modismog('contrast','white')

   %%
thresh =50;

ind = err>thresh;

if any(ind)

   plot([A(ind).X],[A(ind).Y],'r.','markersize',1)
  % title(num2str(kd))
%    axis image
%    axis([axis + 5e3*[-1 1 -1 1]])
%    modismog('contrast','white')
else
   disp bad
end



% This script combines all shapefile data from AutoTerm into a single .mat
% file. This script takes about 8 minutes to run on my laptop. 
% Chad A. Greene, Nov 2022. 

%d = dir('/Users/cgreene/Documents/data/coastlines/AutoTerm/*.shp');


load('autoterm_error_thresholds.mat')

clear A ta
k=1; 

for kd = 1:295
   S_tmp = m_shaperead(['/Users/cgreene/Documents/data/coastlines/AutoTerm/',autoterm_filename{kd}]);
   t_tmp = S_tmp.dbfdata(:,1);
   s_tmp = S_tmp.dbfdata(:,3);
   e_tmp = S_tmp.dbfdata(:,4);


   for sk = 1:length(t_tmp)
   
      tmp = str2double(e_tmp{sk});
      if tmp<error_threshold(kd)

         if kd==124 % manual intervention for alison
            if tmp<59
               continue
            end
         end

         if kd==190 % manual intervention for hayes N'
            if tmp<50
               continue
            end
         end

         if kd==202 % manual intervention 
            if tmp<12.5
               continue
            end
            if tmp>14 & tmp<15
               continue
            end
         end

         if kd==199 % manual intervention for jakobshavn
            if tmp>69 & tmp<70
               continue
            end
            if tmp>124 & tmp<125
               continue
            end
         end

         if kd==212 % manual intervention for kjer
            if tmp<100
               continue
            end
         end

         if kd==279 % manual intervention for dodge
            if tmp<21.3
               continue
            end
         end

         [A(k,1).X,A(k,1).Y] = ll2psn(S_tmp.ncst{sk}(:,2),S_tmp.ncst{sk}(:,1)); 
      
         ta(k,1) = datenum(t_tmp{sk});
         sensor{k,1} = s_tmp{sk}; 
         err(k,1) = tmp; 
         k=k+1; 
      end
   end
   kd
end

for k=1:length(A) 
   A(k).X = A(k).X'; 
   A(k).Y = A(k).Y'; 
end

readme = 'This is AutoTerm data (Zhang 2022, https://zenodo.org/record/7190740) repackaged from 295 shapefiles and trimmed to high-confidence picks to a single .mat file by autoterm_data_reformatter.m.';

% save('/Users/cgreene/Documents/GitHub/greenland-coastlines/data/autoterm_reformatted_trimmed.mat','A','ta','sensor','err','readme','-v7.3')

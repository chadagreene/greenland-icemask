function land = remove_icebergs(land,marine,force_land,force_ocean)
% remove_icebergs finds bits of ice that detach from the ice sheet, and
% removes them from a land mask cube. 
% 
%% Syntax 
% 
%  landr = remove_icebergs(land,marine)
%  landr = remove_icebergs(land,marine,force_land,force_ocean)
% 
%% Description 
% 
% landr = remove_icebergs(land,marine) removes from the logical 3D data cube 
% land any bits of ice that detach from the ice sheet at some point in time
% (assuming dimension 3 of land corresponds to time). The icebergs must be
% within the marine mask (bed<0) to be removed. 
% 
% landr = remove_icebergs(land,marine,force_land,force_ocean) tries to
% ensure that any force_land pixels (a 2D logical mask) are always true,
% and any force_ocean pixels are always false. 

%% Author Info 
% Chad A. Greene, NASA Jet Propulsion Laboratory, August 2022. 

%% Input checks: 

narginchk(2,4)
assert(islogical(land),'Input land must be logical')
assert(islogical(marine),'Input marine must be logical.') 
assert(ndims(land)==3,'Input land must be a 3D cube.') 
assert(isequal([size(land,1),size(land,2)],size(marine)),'Dimensions of land must match dimensions of marine.') 
if nargin>2
   assert(islogical(force_land) & islogical(force_ocean),'Inputs force_land and force_ocean must be logical.')
   assert(isequal(size(force_land),size(force_ocean),size(marine)),'Inputs force_land and force_ocean must be the same size as marine.')
end

%%

[landmasses,n] = bwlabel(all(land,3) & ~marine);

% Get a row,column for each persistent landmass 
r = nan(n,1);
c = r; 
for lm = 1:n
   [r(lm),c(lm)] = find(landmasses==lm,1,'first'); 
end

% Force land or ocean pixels if requested: 
if nargin>=4
   for k = 1:size(land,3)
      tmp = land(:,:,k); 
      tmp(force_land) = true; 
      tmp(force_ocean) = false; 
      land(:,:,k) = tmp;
   end
end

% Select only the pixels that are connected to persistent landmasses: 
for k = 1:size(land,3)
   land(:,:,k) = imfill(bwselect(land(:,:,k),c,r),4,'holes'); 
end
   
% One final time:
if nargin>=4
   for k = 1:size(land,3)
      tmp = land(:,:,k); 
      tmp(force_land) = true; 
      tmp(force_ocean) = false; 
      land(:,:,k) = tmp;
   end
end

end
function [land,landsum] = terminus_masker_advect(x_term,y_term,t_term,x,y,t,land,rock,vx,vy,t_reach,streamreach)
% terminus_masker_advect adjusts a binary ice mask cube to match terminus
% position observations, and fills in un-observed times using known
% velocities. This function operates by artifically advecting terminus
% positions upstream as we go backward in time, declaring that everything 
% upstream of the advected positions must be land (or ice). Then, after
% working backwards, the function repeats itself in a forward manner,
% advecting terminus positions downstream and forward in time, declaring
% that everything downstream of the advected positions must be ocean. 
% 
%% Syntax 
% 
%  land_adj = terminus_masker_advect(x_term,y_term,t_term,x,y,t,land,rock,vx,vy,t_reach)
%  land_adj = terminus_masker_advect(x_term,y_term,t_term,x,y,t,land,rock,vx,vy,t_reach,streamreach)
%  [land_adj,landsum] = terminus_masker_advect(...)
% 
%% Inputs 
% 
% x_term, y_term, t_term: Scattered terminus positions and times (datenum)
% after densification by terminus_data_densifier.m. 
% 
% x, y: Vectors of grid cell center coordinates corresponding to land, vx,
% and vy. 
% 
% t: Datenum times (probably monthly) corresponding to dimension 3 of the
% land cube. 
% 
% land: Binary cube containing true for all ice or rock, false otherwise. 
% 
% rock: 2D binary mask containing true for all rock. 
% 
% vx, vy: 2D velocity grid vectors corresponding to grid cell centers x,y. 
% 
% t_reach: Scalar or vector containing number of days to look ahead and behind
% for terminus positions. For example if t_reach = [500 100], the function
% starts with the last value of t and finds all x_term, y_term points whose
% t_term is less than 500 days after t(end). The points are advected
% backward in time to t(end) and masking is performed. The process is
% repeated for t(end-1) and so on. After getting to t(1), everything is
% repeated in a forward manner, finding the points less than 500 days in
% the past from t(1), advecting them up to expected positions at t(1), then
% masking out all the points downstream of the advected positions. After
% getting to t(end) with t_reach=500, everything is repeated for t_reach=100. 
% 
% streamreach: Length of streamlines in meters. Longer streamlines ensure
% nothing gets missed, but they also take longer to solve and require more
% memory. If streamreach is not specified, the default value is 30e3. 
% 
%% Outputs: 
% 
% land_adj: land cube, adjusted to match terminus position observations. 
% 
% landsum: Optional output provides a time series of the sum of land pixels 
% after each t_reach step. The landsum output is mainly for validating code. 
% 
%% Author Info
% Written by Chad A. Greene of NASA Jet Propulsion Laboratory. 
% July 2022. 

%% Error checks

narginchk(11,12)
assert(isequal(size(x_term),size(y_term),size(t_term)),'Dimensions of x_term, y_term, and t_term must all agree.')
assert(isvector(x) & isvector(y) & isscalar(unique(diff(x))) & isscalar(unique(diff(y))),'Coordinate vectors x and y must be equally spaced 1d arrays.') 
assert(isequal(size(land),[numel(y) numel(x) numel(t)]),'Dimensions of land cube must correspond to lengths of y, x, and t, respectively.')
assert(isequal(size(vx),size(vy),size(rock)),'Dimensions of vx,vy and rock must all be equal')
assert(size(vx,1)==numel(y) & size(vx,2)==numel(x),'Dimensions of vx,vy must match x and y arrays.')
assert(islogical(land),'Input land cube must be binary.')
assert(islogical(rock),'Input rock mask must be binary.')
assert(isvector(t_reach),'t_reach must be a scalar or a vector.') 
assert(~any(diff(t)<=0),'Values of t should be monotonically increasing.')
assert(~any(diff(t_reach)>0),'Values of t_reach should go from big to small.')

%% Enter preferences: 
% Streamline properties define how far upstream or downstream from the
% terminus the streamlines are calculated. The step size is how many steps
% per grid cell are calculated. We want a step smaller than 0.5 to
% essentially meet Nyquist, and a super high density (e.g, step of 0.0001) 
% would ensure that at even a streamline passing through a tiny corner of a
% grid cell will be accurately represented in that grid cell. However,
% small step sizes take longer to solve, so there's a tradeoff. 
% 
% The reach sets a limit on how far each streamline goes. A long reach
% ensures nothing is missed, but double the reach means double the memory
% and/or processing time. 

step = 0.2; % Number of streamline steps per grid cell.  
if nargin<12
   streamreach = 30e3; % meters of reach upstream or downstream from measured termini
end

%% Perform calculations

Npts = round(streamreach/(diff(x(1:2))*step)); % Number of streamline points in each direction 

% Preallocate landsum: 
if nargout==2
   landsum = nan(numel(t),numel(t_reach)); 
end

%disp(['t_range 0 of ',num2str(length(t_reach)),'. ',datestr(now)])

% Do long reaches, then shorter and shorter: 
for kt = 1:length(t_reach)

   % Starting positions of terminus data: 
   x_adv = x_term; 
   y_adv = y_term; 
   t_adv = t_term; 

   % Start with the last timestep and look for future terminus positions: 
   for k = length(t):-1:1

      % Indices of terminus points collected less than t_reach(kt) days into the future, but not in the past:    
      ind = t_adv>=t(k) & t_term<=(t(k)+t_reach(kt)); 

      % Advect the future points upstream to where they should be today:
      dx = interp2(x,y,vx,x_adv(ind),y_adv(ind)) .* (t_adv(ind)-t(k))/365.25;
      dy = interp2(x,y,vy,x_adv(ind),y_adv(ind)) .* (t_adv(ind)-t(k))/365.25;
      x_adv(ind) = x_adv(ind) - dx; 
      y_adv(ind) = y_adv(ind) - dy; 
      t_adv(ind) = t(k); % the points have now been "advected" to "today". 

      xtmp = x_adv(ind); 
      ytmp = y_adv(ind); 
      
      % Ensure no points have advected outside the domain:  
      good = xtmp>x(1) & xtmp<x(end) & ytmp>min(y) & ytmp<max(y); 
      xtmp = xtmp(good); 
      ytmp = ytmp(good); 

      % Figure out which points are already ice, because we won't need to calculate their 
      % backwards streamlines (because that would just end up setting "ice" to "ice"): 
      icy = interp2(x,y,land(:,:,k),xtmp,ytmp,'nearest'); 

      % Upstream streamline calculation: 
      st = stream2(x,y,-vx,-vy,xtmp(~icy),ytmp(~icy),[step Npts]);
      us = cell2mat(st(:)); % stream paths upstream of the terminus
      if ~isempty(us)
         
         % Figure out which streamline points are finite: 
         isf = isfinite(us(:,1)); 
         
         % Grid up the finite streamline points: 
         Us = gridbin(us(isf,1),us(isf,2),true(size(us(isf,1))),x,y,@any); % gridded cells upstream of the terminus 

         % Make a temporary land mask to work with, starting with the land mask we already have for this month: 
         tmp = land(:,:,k); 

         % Make adjustments to the mask: 
         tmp(imerode(imfill(imdilate(Us,true(3)),4,'holes'),true(3))) = true; % set upstream points to true (land)
         
         % Any potential holes in the middle of the land should not be considered ocean, so fill them in:  
         tmp = imfill(tmp,4,'holes'); 
         
         % Remove any skinny ice tongues: 
         tmp = bwmorph(tmp,'open'); 

         % Redefine today's land mask as our adjusted mask: 
         land(:,:,k) = tmp; 
      end

   end

   % REPEAT THE ABOVE, BUT MOVING FORWARD THIS TIME: 
   
   % Reset advection data: 
   x_adv = x_term; 
   y_adv = y_term; 
   t_adv = t_term; 

   % Move forward in time: 
   for k = 1:length(t)

      % Indices of terminus points collected less than t_reach(kt) days ago:    
      ind = t_adv<=t(k) & t_term>=(t(k)-t_reach(kt)); 

      % Advect the previous points downstream to where they should be today:
      dx = interp2(x,y,vx,x_adv(ind),y_adv(ind)) .* (t_adv(ind)-t(k))/365.25;
      dy = interp2(x,y,vy,x_adv(ind),y_adv(ind)) .* (t_adv(ind)-t(k))/365.25;
      x_adv(ind) = x_adv(ind) - dx; 
      y_adv(ind) = y_adv(ind) - dy; 
      t_adv(ind) = t(k); % the points have now been advected to "today". 

      xtmp = x_adv(ind); 
      ytmp = y_adv(ind); 
      
      % Ensure no points have advected outside the domain:  
      good = xtmp>x(1) & xtmp<x(end) & ytmp>min(y) & ytmp<max(y); 
      xtmp = xtmp(good); 
      ytmp = ytmp(good); 

      % If it's already ocean, then the forward streamlines won't affect anything, so don't waste resources calculating their streamlines.  
      icy = interp2(x,y,land(:,:,k),xtmp,ytmp,'nearest'); 

      % Downstream: 
      st = stream2(x,y,vx,vy,xtmp(icy),ytmp(icy),[step Npts]);
      ds = cell2mat(st(:)); % stream paths upstream of the terminus
      if ~isempty(ds)
         
         % Figure out which streamline points are finite: 
         isf = isfinite(ds(:,1)); 
         
         % Grid up the finite streamline points: 
         Ds = gridbin(ds(isf,1),ds(isf,2),true(size(ds(isf,1))),x,y,@any); % gridded cells downstream of the terminus 
         
         % Make a temporary land mask to work with, starting with the land mask we already have for this month: 
         tmp = land(:,:,k); 

         % Make adjustments to the mask: 
         tmp(imerode(imfill(imdilate(Ds,true(3)),4,'holes'),true(3)) & ~rock) = false; % set downstream points to false (ocean)
         
         tmp(rock) = true; % This should already be true, but just in case.
         
         % Any potential holes in the middle of the land should not be considered ocean, so fill them in: 
         tmp = imfill(tmp,4,'holes'); 

         % Remove any skinny ice tongues: 
         tmp = bwmorph(tmp,'open'); 
         
         % Redefine today's land mask as our adjusted mask: 
         land(:,:,k) = tmp; 

      end

   end
   
   if nargout==2
      landsum(:,kt) = squeeze(sum(sum(land))); 
   end
   
   %disp(['t_range ',num2str(kt),' of ',num2str(length(t_reach)),'. ',datestr(now)])
end

%disp(['finished ',datestr(now)])

end % end main function


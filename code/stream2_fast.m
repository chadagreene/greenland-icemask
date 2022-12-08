function [xout,yout] = stream2_fast(x,yf,vxf,vyf,startx,starty,options)

assert(yf(end)>=yf(1),'Flip the order of y.') 
res = x(2)-x(1); 

xind = (startx -  x(1))/res + 1;
yind = (starty - yf(1))/res + 1;

step = options(1); 
maxvert = options(2); 
N = numel(startx); 
xy = cell(1,N); 
for k = 1:N
   xy{k} = stream3c([],[],[],vxf,vyf,[],xind(k),yind(k),[],step,maxvert)';
end

ds = cell2mat(xy(:)); % stream paths from terminus 

xout = (ds(:,1)-1)*res + x(1); 
yout = (ds(:,2)-1)*res + yf(1);

% if ~isempty(ds)
%    isf = isfinite(ds(:,1)); 
%    Ds = gridbin(ds(isf,1),ds(isf,2),true(size(ds(isf,1))),x,y,@any); % Grid cells that are downstream of the terminus 
% else
%    Ds = false(size(X)); 
% end


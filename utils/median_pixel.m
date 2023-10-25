function medval = median_pixel(x,y,val,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% medval = median_pixel(x,y,val)
%
% Imagine you have a set of data with 'val' that is distributed at locations
% x,y. Given that a unique location (x,y) might have several data points, you
% want to have a stable estimate (median/mean) of the value of each unique
% location. 
% --output 'MEDVAL' now has 7 columns: [x y num median/mean min max std], so the
%   the first 3 columns are the same as in 'density_pixel'. Use 'density_pixel'
%   if you are interested in counts a unique locations only.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/10/05
% Last modified date:   2022/10/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('method','median');

x = reshape(x,[],1);  % force to be a column vector
y = reshape(y,[],1);
val = reshape(val,[],1);
loc = [x y]; 

%find the unique entries
locunq = unique(loc,'rows','stable');
medval = zeros(size(locunq,1), 7);
for i = 1: size(locunq,1)
  ind = find(ismember(loc ,locunq(i,:), 'rows'));
  medval(i,1:2) = locunq(i,:);
  if isequal(method, 'median')
    medval(i,3) = median(val(ind));
  elseif isequal(method, 'mean')
    medval(i,3) = mean(val(ind));
  end
  medval(i,4) = length(ind);
  medval(i,5) = min(val(ind)); 
  medval(i,6) = max(val(ind)); 
  medval(i,7) = std(val(ind));
end

% keyboard
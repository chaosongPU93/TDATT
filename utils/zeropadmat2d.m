function [xyzgridpad,xgrid,ygrid,zgrid,indgrid] = ...
  zeropadmat2d(xyzscatter,xvec,yvec,const)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xyzgridpad,xgrid,ygrid,zgrid,ind2] = zeropadmat2d(xyzscatter,xvec,yvec)
%
% Imagine you have a sparse scattered 2D data set (which has an intrisic 
% sampling rate) with the x and y as locations and z as the some value on
% that location. So the 'xyzscatter' has 3 columns [x y z]. 
% What you want to do now is to zero-padding this data set to a 2D matrix
% where the rectangular grid is defined by the 'xvec' and 'yvec'. Essentially
% the grid points [x,y] that were defined in 'xyzscatter' has the same z 
% value, while ones that were missing are filled with zeros, such that the 
% return 'zgrid' is a 2D matrix. The created grid 'xgrid' and 'ygrid' from 
% 'xvec' and 'yvec'are also returned. The main output 'xyzgridpad' has the
% same format as the input 'xyzscatter' that has 3 colomns.
%
% --2022/10/05, modify from padding 0s to padding any constant as desired, but
%   the default is still zero.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/04/13
% Last modified date:   2022/10/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('const',0);

%force the input grid range to a column vector
xvec = reshape(xvec,[],1);
yvec = reshape(yvec,[],1);
if yvec(1)<yvec(end)
  yvec = flipud(yvec);  % flip the order to make the grid have right-hand positive axes
end

%create a rectangular grid
[xgrid, ygrid] = meshgrid(xvec,yvec);
zgrid = zeros(size(xgrid))+const*ones(size(xgrid));

%make them a column vector for easy processing
x = reshape(xgrid,[],1);
y = reshape(ygrid,[],1);
z = reshape(zgrid,[],1);

%find the index in created grid that correspond to the input scattered x y z
xyloc = [x y];
% [~,ind1,ind2] = intersect(xyzscatter(:,1:2),xyloc,'rows','stable');
% %if not all input has a match of grid point
% if ~isequaln(ind1,(1:size(xyzscatter,1))') || ~isequal(length(ind2),size(xyzscatter,1))   
%   error('The input range of grid is smaller than data scatter.');
% end
ind2 = nan(size(xyzscatter,1),1);
for i = 1:size(xyzscatter,1)
  ind2(i) = find(abs(xyloc(:,1)-xyzscatter(i,1))<1e-6 & abs(xyloc(:,2)-xyzscatter(i,2))<1e-6);
  if isempty(ind2(i))
    error('The input range of grid is smaller than data scatter.');
  end
end

%At those matched grid points, just replace the z value 
z(ind2) = xyzscatter(:,3);

%revert z vector to the grid matrix
zgrid = reshape(z,size(xgrid));

%return the result as the same format as input 'xyzscatter'
xyzgridpad = [x y z];

%return indices of zero-padded grid that correspond to input 'xyzscatter'
indgrid = ind2;

% keyboard












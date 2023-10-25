function [density,induniq] = density_pixel(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [density,indices] = density_pixel(x,y)
%
% This function is to obtain the density (number of points;
% repetition times) of a 2-D data set (x,y). Different from
% 'density_matrix' which bins data based on a bin of 
% rectangle grid, this function bins data pixel wise. It 
% simply counts the times that the same (x,y) repeats. 
% Therefore, the input data set should have identical entries.
% The output 'density' has 3 columns, first 2 of which are 
% x, y locations while the 3rd columns is the counts at this
% location. Also return the indices of data that being binned
% to the same pixel. 
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/14
% Last modified date:   2023/05/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(x,[],1);  % force to be a column vector
y = reshape(y,[],1);
data = [x y]; 

%find the unique entries
dataunq = unique(data,'rows','stable');
nuniq = size(dataunq,1);
density = zeros(nuniq, 3);
induniq = cell(nuniq,1);
for i = 1: nuniq
  ind = find(ismember(data ,dataunq(i,:), 'rows'));
  density(i,1:2) = dataunq(i,:);
  density(i,3) = length(ind);
  induniq{i} = ind;
end





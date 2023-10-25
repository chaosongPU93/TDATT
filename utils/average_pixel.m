function average = average_pixel(x,y,valmat,avemethod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average = average_pixel(x,y,valmat,avemethod)
%
% This function is to obtain the average of an property come with a data set
% at each location x, y. In other words, 'valmat' is the value of some property
% at each disrecte point x, y. Could have multiple properties as different 
% columns.
% We want to know what is the average value of it
% at each unique location. Can choose one of 2 methods of averaging, mean or
% median. Similar to 'density_pixel' which just counts the number of unique
% locations
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/07/21
% Last modified date:   2022/07/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(x,[],1);  % force to be a column vector
y = reshape(y,[],1);
data = [x y]; 
nval = size(valmat,2);   % how many types of properties?

%find the unique entries
dataunq = unique(data,'rows','stable');
average = zeros(size(dataunq,1), 2+nval);
for i = 1: size(dataunq,1)
  average(i,1:2) = dataunq(i,:);
  if isequal(avemethod, 'mean') 
    average(i,3:end) = mean(valmat(ismember(data ,dataunq(i,:), 'rows'), :), 1);
  elseif isequal(avemethod, 'median') 
    average(i,3:end) = median(valmat(ismember(data ,dataunq(i,:), 'rows'), :), 1);
  end
end
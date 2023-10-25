function [xcnt,ycnt,x1sig,y1sig] = ranybinxeqnum(x,y,avgmethod,nbin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xcnt,ycnt,x1sig,y1sig] = ranybinxeqnum(x,y,avgmethod,nbin)
%
% This function tries to do something similar to 'ranybinx.m', but instead
% of setting bins of equal binwidth, it sets bins where the number of data
% that falls into each bin is the same. Ideally, if each member of the data
% set has the same 'error', then the error bar of each bin with the same 
% number of members should be the same as well, statistically. 
% 
% --Note that you can only specify 'nbin'. The data sorted in x dimension 
%   then would be divided into nbins, and a median or mean based on the
%   averaging method would be obtained for all members in that bin.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/07/05
% Last modified date:   2023/07/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('avgmethod','median');
% defval('nbin',10);

%force to a column vector
x = reshape(x,[],1);  
y = reshape(y,[],1);
ntot = length(x);
n = round(ntot/nbin);

%sort data
[xsort,isort] = sort(x,'ascend');
ysort = y(isort);

xcnt = zeros(nbin,1);
ycnt = zeros(nbin,1);
x1sig = zeros(nbin,1);
y1sig = zeros(nbin,1);
for i = 1: nbin
  if i~= nbin
    xbin = xsort((i-1)*n+1: i*n);
    ybin = ysort((i-1)*n+1: i*n);
  else
    xbin = xsort((i-1)*n+1: end);
    ybin = ysort((i-1)*n+1: end);
  end
  if isequal(avgmethod, 'median')
    xcnt(i) = median(xbin);
    ycnt(i) = median(ybin);
  elseif isequal(avgmethod, 'mean')
    xcnt(i) = mean(xbin);
    ycnt(i) = mean(ybin);
  end
  x1sig(i) = std(xbin);
  y1sig(i) = std(ybin);
end
  
% keyboard
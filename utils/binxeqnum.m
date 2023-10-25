function [xbin,indbin,nxbin] = binxeqnum(x,nbin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xbin,indbin,nxbin] = binxeqnum(x,nbin)
%
% This function sets 'nbin' bins where the number of data 'x'
% that falls into each bin is the same.
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

defval('nbin',10);

%force to a column vector
x = reshape(x,[],1);  
ntot = length(x); %total num
n = round(ntot/nbin); %num in each bin 

%sort data
[xsort,indsort] = sort(x,'ascend');

xbin = cell(nbin,1);  %x data in each bin
indbin = cell(nbin,1);  %indices of x data in each bin
nxbin = zeros(nbin,1);  %actual num of x data in each bin
for i = 1: nbin
  if i~= nbin
    xbin{i} = xsort((i-1)*n+1: i*n);
    indbin{i} = indsort((i-1)*n+1: i*n);
  else
    xbin{i} = xsort((i-1)*n+1: end);
    indbin{i} = indsort((i-1)*n+1: end);
  end
  nxbin(i) = length(xbin{i});
end
  
% keyboard
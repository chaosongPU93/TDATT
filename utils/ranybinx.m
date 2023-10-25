function [xcnt,ycnt,y1sig] = ranybinx(x,y,avgmethod,nbin,binw,binEdges)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xcnt,ycnt,y1sig] = ranybinx(x,y,avgmethod,nbin,binw)
%
% Imagine you have a dataset which has properties 'x' and 'y', between which
% you want to know if there is any relation. The simplest way is just to plot
% the scatter x and y. This function does one more step. It tries to bin the 
% property 'x' based on the number of bins 'nbin', OR, bin width 'binw'. For
% the data points whose property 'x' fall into that bin, obtain the median or 
% mean of property 'y', and its +-1\sigma range. 
% 
% --Note that you should specify either 'nbin' or 'binw', but not both of them,
%   the other one is computed accordingly.
% --previously bin center/edge is determined by actual range of x. Now you can
%   specify binedge directly so that different data sets can be compared based
%   on common bins.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/28
% Last modified date:   2022/11/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('avgmethod','median');
% defval('nbin',10);
% defval('binw',[]);
% defval('binEdges',[]);

%force to a column vector
x = reshape(x,[],1);  
y = reshape(y,[],1);

%if 'binw' is specified
if isempty(nbin) && ~isempty(binw) && isempty(binEdges)
  nbin = ceil((max(x)-min(x))/binw);
  binEdges = min(x): binw: min(x)+nbin*binw;
%if 'nbin' is specified 
elseif ~isempty(nbin) && isempty(binw) && isempty(binEdges)
  binw = (max(x)-min(x))/nbin;
  binEdges = min(x): binw: min(x)+nbin*binw;
%if 'binEdges' is specified 
elseif isempty(nbin) && isempty(binw) && ~isempty(binEdges)
  binw = binEdges(2)-binEdges(1);
  nbin = length(binEdges)-1;
%otherwise return error
else
  error("You can only specify 'nbin' or 'binw' or 'binEdges' \n");
end
    
xcnt = zeros(nbin,1);
ycnt = zeros(nbin,1);
y1sig = zeros(nbin,1);
for i = 1: nbin
  ybin = y(x>=binEdges(i) & x<binEdges(i+1)); % y fall into the bin
  xcnt(i) = (binEdges(i)+binEdges(i+1))/2;
  if ~isempty(ybin)
    if isequal(avgmethod, 'median')
      ycnt(i) = median(ybin);
    elseif isequal(avgmethod, 'mean')
      ycnt(i) = mean(ybin);
    end
    y1sig(i) = std(ybin);  
  else
    ycnt(i) = nan;
    y1sig(i) = nan;
  end
end

% keyboard
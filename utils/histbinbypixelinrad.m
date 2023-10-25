function [bincenter,binheight,count,normalizer] = ...
  histbinbypixelinrad(datavect,binedge,normalization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = axsym(ax,w)
%
% This function is to bin the 'datavect' by 2D circular area in terms of 
% integer number of pixel inside the area circumvented with a radius, defined
% by a series of concentric circles. Different from 'histbinbyarea' which 
% computes the actual area, this just count the integer pixels that fall in 
% the bin edges, which will be used as the normalizer for the counts. 
% Note that the convention is that left side of the bin edge is always closed.
% For example, if the bin edge is like 0, 1, 2, etc. Then 0 is 
% put into bin 0-1, 1 is put into bin 1-2.
% 
% 
%         
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/21
% Last modified date:   2023/03/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bincenter = (binedge(2:end) + binedge(1:end-1))/2;
nbin = length(binedge)-1;
count = zeros(nbin,1);
npix = zeros(nbin,1);
for i = 1: nbin
  count(i) = sum(datavect>=binedge(i) & datavect<binedge(i+1));
  for ii = -binedge(i+1): 1: binedge(i+1)
    for jj = -binedge(i+1): 1: binedge(i+1)
      if sqrt(ii^2+jj^2)<binedge(i+1)
        npix(i) = npix(i)+1;
      end
    end
  end
end  
npixdiff = [npix(1); diff(npix)];

if strcmp(normalization,'count')
  normalizer = ones(nbin,1);
elseif strcmp(normalization,'countdensity')
  normalizer = npixdiff; % / binedge(end)^2)
elseif strcmp(normalization,'probability')
  normalizer = length(datavect) * ones(nbin,1);
elseif strcmp(normalization,'pdf')
  normalizer = length(datavect) * npixdiff;
end

binheight = count./normalizer;
  
% keyboard  
  
  
  
  
  
  
  
  
  
  
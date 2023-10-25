function [bincenter,binheight,count,normalizer] = ...
  histbinbyarea(datavect,binedge,normalization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [bincenter,binheight,count,normalizer] = ...
%   histbinbyarea(datavect,binedge,normalization)
%
% This function is to bin the 'datavect' by 2D circular area. Imagine data
% is like some 2D Euclidean distance, due to its sqrt root nature, 1D histogram
% along each dimension separately could not best visualize its distribution.
% Hence here we can try to bin it by the circular area. Using a series of 
% radius with a fixed increment 'dr', we can define 'nbin' of area between
% two consecutive circles. 
% Then for bin i, where i=1,2,..., area Ai = \pi * [ (i*dr)^2 - ((i-1)*dr)^2 ]
% Note that the convention is that left side of the bin edge is always closed.
% For example, if the bin edge is like 0, 1, 2, etc. Then 0 is 
% put into bin 0-1, 1 is put into bin 1-2.
% 
%         
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/03/13
% Last modified date:   2023/03/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% binedge = (0: dr: nbin*dr)';
binedge = reshape(binedge,[],1);
bincenter = (binedge(2:end) + binedge(1:end-1))/2;
nbin = length(binedge)-1;
count = zeros(nbin,1);
% normalizer = zeros(nbin,1);
% binheight = zeros(nbin,1);
for i = 1: nbin
  count(i) = sum(datavect>=binedge(i) & datavect<binedge(i+1));
end

if strcmp(normalization,'count')
  normalizer = ones(nbin,1);
elseif strcmp(normalization,'countdensity')
  normalizer = pi*(binedge(2:end).^2 - binedge(1:end-1).^2); % / binedge(end)^2)
elseif strcmp(normalization,'probability')
  normalizer = length(datavect) * ones(nbin,1);
elseif strcmp(normalization,'pdf')
  normalizer = length(datavect) * pi*(binedge(2:end).^2 - binedge(1:end-1).^2);
end

binheight = count./normalizer;


% keyboard
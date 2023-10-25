function [newt,trange,ntot] = burst_range(burst,intertime,refjuldate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to convert the index of the detections inside each tremor
% burst to the occurrence times in year, julday and sec, the same format
% as 'trange'.
%
% newt: is simply the start and end time of each burst, in relative time of
%       days
% trange: is times in year+julday and the start and end time of each burst
%         in sec, the same format as 'trange'
% ntot: total number of detections in all bursts
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/26
% Last modified date:   2022/03/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newt = zeros(size(burst,1), 2);
trange = zeros(size(burst,1), 3);
ntot = 0;
for i = 1: size(burst,1)
    ntot = ntot+size(burst{i},1);
    ind = burst{i};
    newt(i,1) = intertime(ind(1),1); %-0.5/60/24
    newt(i,2) = intertime(ind(end),1); %+0.5/60/24
end
trange(:,1) = floor(newt(:,1))+refjuldate;
trange(:,2:3) = (newt(:,1:2)-floor(newt(:,1))).*(3600.*24);

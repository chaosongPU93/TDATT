function [septime] = interevt_time(occtime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [septime] = interevt_time(occtime)
% This function is to read in the time of detections and return the 
% separation in time between itself and its preceding detection,
% i.e., inter-event time
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/26
% Last modified date:   2021/09/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

septime = zeros(length(occtime),2);
% 1st col: occurence time
% 2nd col: inter-detection time to its preceding detection
septime(:,1) = occtime;
for i = 1: size(septime,1)
    if i == 1
%         septime(i,2) = septime(i+1,1)-septime(i,1); % 1st detection has no preceding, so use time to its following  
        septime(i,2) = 0; % 1st detection has no preceding, so use 0
    else
        septime(i,2) = septime(i,1)-septime(i-1,1);
    end
end

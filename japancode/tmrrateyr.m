function [rateyr] = tmrrateyr(yrall,tmrall,narea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to calculate the tremor detection rate every year within 
% the input year range
%   
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/24
% Last modified date:   2020/02/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rateyr = [];

for i = 1: length(yrall)
% i=length(yrall);
    if i == 1
        date1 = tmrall(1,1);
    else
        date1 = yrall(i)*10000 + 101;
    end
    
    if i == length(yrall)
        date2 = tmrall(end,1);
    else
        date2 = yrall(i)*10000 + 1231;
    end
        
    d1 = datetime(date1,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    d2 = datetime(date2,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    
    % number of days in total in that year
    ndyyr(i) = caldays(between(d1,d2,'days'))+1;
    
    % number of events in toal in that year
    tmrallyr = tmrall(tmrall(:,1)>=date1 & tmrall(:,1)<=date2, :);
    nallyr(i) = size(tmrallyr,1);
    
    % occurence rate, num per day per unit area
    rateyr(i) = nallyr(i)/ndyyr(i)/narea;
    
end


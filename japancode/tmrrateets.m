function [rttmrets] = tmrrateets(stday,edday,tmrall,narea)
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

d1  = datetime(tmrall(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

for i = 1: length(stday)
    % start and end date of ETS 
    dis = datetime(d1+caldays(stday(i))-1,'Format','yyyy-MM-dd');
    die = datetime(d1+caldays(edday(i))-1,'Format','yyyy-MM-dd');
    dateis = yyyymmdd(dis');
    dateie = yyyymmdd(die');
    
    ndyets(i) = edday(i)-stday(i)+1;
    % tremor detections during ETS
    tmrets = tmrall(tmrall(:,1)>=dateis & tmrall(:,1)<=dateie, :);
    ntmrets(i) = size(tmrets,1);
    rttmrets(i) = ntmrets(i)/ndyets(i)/narea;
    
end

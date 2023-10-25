function rateyrm = evtrateyr(yrall,evtall,tmrall,narea,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to calculate the event occurence rate per day per unit area
% every year within the input year range, and classify further according to 
% the magnitude  
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/24
% Last modified date:   2020/02/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rateyr = [];
rateyrm0 = [];
rateyrm01 = [];
rateyrm12 = [];
rateyrm23 = [];
rateyrm34 = [];
rateyrm4 = [];

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
    ndayyr(i) = caldays(between(d1,d2,'days'))+1;
    
    % number of events in toal in that year
    evtallyr = evtall(evtall(:,1)>=date1 & evtall(:,1)<=date2, :);
    % all magnitude events
    nallyr(i) = size(evtallyr,1);
    rateyr(i) = nallyr(i)/ndayyr(i)/narea;
    
    if class == 7
        % events with different magnitudes
        nallyrm0(i) = size(evtallyr(evtallyr(:,11)<0, :), 1);
        nallyrm01(i) = size(evtallyr(evtallyr(:,11)>=0 & evtallyr(:,11)<1, :), 1);
        nallyrm12(i) = size(evtallyr(evtallyr(:,11)>=1 & evtallyr(:,11)<2, :), 1);
        nallyrm23(i) = size(evtallyr(evtallyr(:,11)>=2 & evtallyr(:,11)<3, :), 1);
        nallyrm34(i) = size(evtallyr(evtallyr(:,11)>=3 & evtallyr(:,11)<4, :), 1);
        nallyrm4(i) = size(evtallyr(evtallyr(:,11)>=4, :), 1);
        % occurence rate, num per day per unit area
        rateyrm0(i) = nallyrm0(i)/ndayyr(i)/narea;
        rateyrm01(i) = nallyrm01(i)/ndayyr(i)/narea;
        rateyrm12(i) = nallyrm12(i)/ndayyr(i)/narea;
        rateyrm23(i) = nallyrm23(i)/ndayyr(i)/narea;
        rateyrm34(i) = nallyrm34(i)/ndayyr(i)/narea;
        rateyrm4(i) = nallyrm4(i)/ndayyr(i)/narea;

    elseif class == 3
        nallyrml1(i) = size(evtallyr(evtallyr(:,11)<1, :), 1);
        nallyrmg1(i) = size(evtallyr(evtallyr(:,11)>=1, :), 1);
        
        rateyrml1(i) = nallyrml1(i)/ndayyr(i)/narea;
        rateyrmg1(i) = nallyrmg1(i)/ndayyr(i)/narea;
    end
    
    
end
if class == 7
    rateyrm = [rateyr;rateyrm0;rateyrm01;rateyrm12;rateyrm23;rateyrm34;rateyrm4];
elseif class == 3
    rateyrm = [rateyr;rateyrml1;rateyrmg1];
end
    
    
    
    
    
    

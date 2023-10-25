function [rttmryrm,rtfreeyrm] = evtrateintmryr(yrall,evtall,tmrall,datetmr,narea,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to calculate the event occurence rate per day per unit area
% in tremor days and tremor free days (depends on the dates fed) in every
% year within the input year range, and classify further according to the
% magnitude 
%   
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/25
% Last modified date:   2020/02/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1: length(yrall)
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
    
    % dates of tremor in each year
    datetmryr = datetmr(datetmr>=date1 & datetmr<=date2);
    % number of days of tremor in each year
    ntmrdyyr(i) = length(datetmryr);
    % number of days in total in that year
    ndyyr(i) = caldays(between(d1,d2,'days'))+1;
    % number of free days  in each year
    nfreedyyr(i) = ndyyr(i)-ntmrdyyr(i);

    % number of events in toal in that year
    evtallyr = evtall(evtall(:,1)>=date1 & evtall(:,1)<=date2, :);
    nallyr(i) = size(evtallyr,1);
       
    % all magnitude events during tremor days 
    evttmryr = [];
    for j = 1: ntmrdyyr(i)
        tmp = evtallyr(evtallyr(:,1)==datetmryr(j), :);
        evttmryr = [evttmryr; tmp];
    end
    ntmryr(i) = size(evttmryr,1);
    rttmryr(i) = ntmryr(i)/ntmrdyyr(i)/narea;
    if class == 7
        % events with different magnitudes
        ntmryrm0(i) = size(evttmryr(evttmryr(:,11)<0, :), 1);
        ntmryrm01(i) = size(evttmryr(evttmryr(:,11)>=0 & evttmryr(:,11)<1, :), 1);
        ntmryrm12(i) = size(evttmryr(evttmryr(:,11)>=1 & evttmryr(:,11)<2, :), 1);
        ntmryrm23(i) = size(evttmryr(evttmryr(:,11)>=2 & evttmryr(:,11)<3, :), 1);
        ntmryrm34(i) = size(evttmryr(evttmryr(:,11)>=3 & evttmryr(:,11)<4, :), 1);
        ntmryrm4(i) = size(evttmryr(evttmryr(:,11)>=4, :), 1);
        % occurence rate, num per day per unit area
        rttmryrm0(i) = ntmryrm0(i)/ntmrdyyr(i)/narea;
        rttmryrm01(i) = ntmryrm01(i)/ntmrdyyr(i)/narea;
        rttmryrm12(i) = ntmryrm12(i)/ntmrdyyr(i)/narea;
        rttmryrm23(i) = ntmryrm23(i)/ntmrdyyr(i)/narea;
        rttmryrm34(i) = ntmryrm34(i)/ntmrdyyr(i)/narea;
        rttmryrm4(i) = ntmryrm4(i)/ntmrdyyr(i)/narea;
        
    elseif class == 3
        ntmryrml1(i) = size(evttmryr(evttmryr(:,11)<1, :), 1);
        ntmryrmg1(i) = size(evttmryr(evttmryr(:,11)>=1, :), 1);
        rttmryrml1(i) = ntmryrml1(i)/ntmrdyyr(i)/narea;
        rttmryrmg1(i) = ntmryrmg1(i)/ntmrdyyr(i)/narea;
    end
    
    % all magnitude events during tremor free days 
    evtfreeyr = setdiff(evtallyr, evttmryr,'rows');
    nfreeyr(i) = size(evtfreeyr,1);
    rtfreeyr(i) = nfreeyr(i)/nfreedyyr(i)/narea;
    if class == 7
        % events with different magnitudes
        nfreeyrm0(i) = size(evtfreeyr(evtfreeyr(:,11)<0, :), 1);
        nfreeyrm01(i) = size(evtfreeyr(evtfreeyr(:,11)>=0 & evtfreeyr(:,11)<1, :), 1);
        nfreeyrm12(i) = size(evtfreeyr(evtfreeyr(:,11)>=1 & evtfreeyr(:,11)<2, :), 1);
        nfreeyrm23(i) = size(evtfreeyr(evtfreeyr(:,11)>=2 & evtfreeyr(:,11)<3, :), 1);
        nfreeyrm34(i) = size(evtfreeyr(evtfreeyr(:,11)>=3 & evtfreeyr(:,11)<4, :), 1);
        nfreeyrm4(i) = size(evtfreeyr(evtfreeyr(:,11)>=4, :), 1);
        % occurence rate, num per day per unit area
        rtfreeyrm0(i) = nfreeyrm0(i)/nfreedyyr(i)/narea;
        rtfreeyrm01(i) = nfreeyrm01(i)/nfreedyyr(i)/narea;
        rtfreeyrm12(i) = nfreeyrm12(i)/nfreedyyr(i)/narea;
        rtfreeyrm23(i) = nfreeyrm23(i)/nfreedyyr(i)/narea;
        rtfreeyrm34(i) = nfreeyrm34(i)/nfreedyyr(i)/narea;
        rtfreeyrm4(i) = nfreeyrm4(i)/nfreedyyr(i)/narea;
    
    elseif class == 3
        nfreeyrml1(i) = size(evtfreeyr(evtfreeyr(:,11)<1, :), 1);
        nfreeyrmg1(i) = size(evtfreeyr(evtfreeyr(:,11)>=1, :), 1);
        rtfreeyrml1(i) = nfreeyrml1(i)/nfreedyyr(i)/narea;
        rtfreeyrmg1(i) = nfreeyrmg1(i)/nfreedyyr(i)/narea;
    end
end

% construct a matrix
if class == 7
    rttmryrm = [rttmryr;rttmryrm0;rttmryrm01;rttmryrm12;rttmryrm23;rttmryrm34;rttmryrm4];
    rtfreeyrm = [rtfreeyr;rtfreeyrm0;rtfreeyrm01;rtfreeyrm12;rtfreeyrm23;rtfreeyrm34;rtfreeyrm4];
elseif class == 3
    rttmryrm = [rttmryr;rttmryrml1;rttmryrmg1];
    rtfreeyrm = [rtfreeyr;rtfreeyrml1;rtfreeyrmg1];
end
    

























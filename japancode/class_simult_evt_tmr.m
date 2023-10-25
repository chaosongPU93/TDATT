function [evtin, evtout] = class_simult_evt_tmr(refdate,tmrobj,evtobj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to plot the difference in minimum distances from the 
% piercing points of events to the closest tremors between 2 velocity models 
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/29
% Last modified date:   2020/04/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% relative time in days for the first day of tremor catalog           
d1 = datetime(refdate,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

% relative origin time in days of all tremors
ditmr = datetime(tmrobj(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ottmr = caldays(between(d1,ditmr,'days')) + tmrobj(:,5)/24;      % unit is day

evtin = [];
evtout = [];
for i = 1:size(evtobj,1)

    % relative origin time in days of all events
    dievt = datetime(evtobj(i,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    otevt = caldays(between(d1,dievt,'days')) + (evtobj(i,5)+evtobj(i,6)/60+evtobj(i,7)/3600)/24;
    
    % travel time difference in day, this difference is theoretically in secs, so is
    % negligible compared with the tremor detection time resolution (hr)
    
    % find all pairs that must pass through the same area at the same time
    otdt = otevt-ottmr;
    if sum(otdt>=0 & otdt<1/24)>=1   % must within the same hr
        evtin = [evtin; evtobj(i,:)];
    else
        evtout = [evtout; evtobj(i,:)];
    end

end


% keyboard
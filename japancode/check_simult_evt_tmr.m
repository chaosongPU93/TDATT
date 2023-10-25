function [indin, indout] = check_simult_evt_tmr(refdate,tmr,evtobj,indtmr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to output the index of tremor that shares the same time as
% the target event (i.e., occurred 1 hr or less before event); and the tremor that
% are outside this time period
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/29
% Last modified date:   2020/07/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

tmrobj = tmr(indtmr,:);

%%% relative time in days for the first day of tremor catalog           
d1 = datetime(refdate,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

% relative origin time in days of all tremors
ditmr = datetime(tmrobj(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ottmr = caldays(between(d1,ditmr,'days')) + tmrobj(:,5)/24;      % unit is day

% relative origin time in days of all events
dievt = datetime(evtobj(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
otevt = caldays(between(d1,dievt,'days')) + (evtobj(:,5)+evtobj(:,6)/60+evtobj(:,7)/3600)/24;

otdt = (otevt-ottmr)*24;   % time difference in hr

indin = indtmr(otdt>=0 & otdt<1);

indout = setdiff(indtmr, indin);



% keyboard
    

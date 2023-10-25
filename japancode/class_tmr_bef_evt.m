function [objtmr, mindt] = class_tmr_bef_evt(refdate,tmr,evtobj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to output the difference in time of the closest tremor
% before the query event 
% i don't need to know the detailed location or other infomation of every
% tremor, all i want is what is cloest tremor in time relative to that query
% event
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/24
% Last modified date:   2020/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

%%% relative time in days for the first day of tremor catalog           
d1 = datetime(refdate,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

% relative origin time in days of all tremors
ditmr = datetime(tmr(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ottmr = caldays(between(d1,ditmr,'days')) + tmr(:,5)/24;      % unit is day

% relative origin time in days of all events
dievt = datetime(evtobj(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
otevt = caldays(between(d1,dievt,'days')) + (evtobj(:,5)+evtobj(:,6)/60+evtobj(:,7)/3600)/24;

otdt = (ottmr-otevt)*24;   % time difference in hr

% tremors that occur after the event are irrelavant, won't affect the data
tmpind = find(otdt < 0);

[mindt, i] = min(abs(otdt(tmpind)));
objind = tmpind(i);
objtmr = tmr(objind);




% keyboard
function [evts,evtti,evtto] = ...
    evtouttmr(evtinall,evtoutall,evtref,tmrall,tmrref,magtol,dloctol,deptol,dttol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a function to find the events that is close to the reference events in
% space within a loction tolerence during the objective dates from a catalog 
%
%
%   INPUT:
%       f:      the construct containing frame and for figure
%       d1:     the start day of the time axis, Must be datestring format
%       d2:     the end day of the time axis, Must be datestring format
%       ievtday:    the regular event day on the time axis, counting from the start 
%       nevtobj:    the event counts on that day
%
%   OUTPUT:
%       f:    the construct for figure
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/20
% Last modified date:   2020/02/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 1. find events that are spatially close to reference, but in tremor free days
% 'timexxx' for datetime type
% 'datexxx' for numerical date like yyyymmdd
% 'datestrxxx' for string type yyyymmdd
% 'idxxx' for the numerical order of one day relative to the first day
% 'ndxxx' for the number of days
time1 = datetime(tmrall(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
time2 = datetime(tmrall(end,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ndall = caldays(between(time1,time2,'days'))+1;

datetmr = unique(tmrall(:,1));
timetmr = datetime(datetmr,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
idtmr = caldays(between(time1,timetmr,'days'))+1;
idfree = setdiff(1:ndall,idtmr);
timefree = datetime(time1+caldays(idfree)-1,'Format','yyyy-MM-dd');
datefree = yyyymmdd(timefree');

% get the events during the objective dates
lend = length(datefree);
evttmp = [];
for i = 1: lend
    evttmp = [evttmp; evtinall(datefree(i) == evtinall(:,1), :)];
end

evts = [];      % events spatially close 
if ~isempty(evttmp)
    for i = 1:size(evtref,1)
        dist = sqrt((evttmp(:,8)-evtref(i,8)).^2 + (evttmp(:,9)-evtref(i,9)).^2);
        evttmp2 = evttmp(dist<=dloctol & abs(evttmp(:,11)-evtref(i,11))<=magtol & ...
                         abs(evttmp(:,10)-evtref(i,10))<=deptol, :);
        % add a marker as the last column to indicate which ref event they correspond to
        evttmp2 = [evttmp2 i*ones(size(evttmp2,1),1)];               
        evts = [evts; evttmp2]; 
    end
end
[~,ievts,~] = unique(evts(:,1:end-1),'rows','stable');
evts = evts(ievts,:);


%%% 2. find event that are timely close to reference, but outside the radius regions

% relative time in days of all reference events
d1 = datetime(tmrall(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
di = datetime(evtref(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
drefrela = caldays(between(d1,di,'days')) + (evtref(:,5)+evtref(:,6)/60+evtref(:,6)/3600)/24;

[indr,~] = boundary(tmrref(:,6),tmrref(:,7),0.5);
[is,ion] = inpolygon(evtinall(:,8),evtinall(:,9),tmrref(indr,6),tmrref(indr,7));
evttmp = evtinall(~is & ~ion, :);

evtti = [];     % events timely close, but outside the tremor region defined by the tmrir
if ~isempty(evttmp)
    % relative time in days of all events outside radius, but still inside the tremor region  
    di = datetime(evttmp(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    dtmprela = caldays(between(d1,di,'days')) + (evttmp(:,5)+evttmp(:,6)/60+evttmp(:,6)/3600)/24;
    for i = 1:size(evtref,1)
        evttmp2 = evttmp(abs(dtmprela-drefrela(i))*24<=dttol & ... 
                         abs(evttmp(:,11)-evtref(i,11))<=magtol & ...
                         abs(evttmp(:,10)-evtref(i,10))<=deptol, :);
        % add a marker as the last column to indicate which ref event they correspond to                     
        evttmp2 = [evttmp2 i*ones(size(evttmp2,1),1)];             
        evtti = [evtti; evttmp2];              
    end
end
[~,ievtti,~] = unique(evtti(:,1:end-1),'rows','stable');
evtti = evtti(ievtti,:);        

%%% 3. find event that are timely close to reference, but outside the tremor regions

evtto = [];     % events timely close, but outside the tremor region defined by the tmrir
if ~isempty(evtoutall)
    % relative time in days of all events outside radius, but still inside the tremor region  
    di = datetime(evtoutall(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    doutrela = caldays(between(d1,di,'days')) + (evtoutall(:,5)+evtoutall(:,6)/60+...
                       evtoutall(:,6)/3600)/24;
    for i = 1:size(evtref,1)
        evttmp2 = evtoutall(abs(doutrela-drefrela(i))*24<=dttol & ... 
                         abs(evtoutall(:,11)-evtref(i,11))<=magtol & ...
                         abs(evtoutall(:,10)-evtref(i,10))<=deptol, :);
        % add a marker as the last column to indicate which ref event they correspond to                     
        evttmp2 = [evttmp2 i*ones(size(evttmp2,1),1)];
        evtto = [evtto; evttmp2];              
    end
end
[~,ievtto,~] = unique(evtto(:,1:end-1),'rows','stable');
evtto = evtto(ievtto,:);
























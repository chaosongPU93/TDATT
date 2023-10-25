% function contemp_evttmr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main script to analyze the spatio-temporal evolution of tremor
% and regular events occurred in the same ETS period, to find the nearly, and/or
% co-located tremor and eqs, save those candidate events for further analysis
% 
%   NOTICE:
%       1. Refer to the record format if not sure about anything
%          for regular events, go to:     
%           https://www.data.jma.go.jp/svd/eqev/data/bulletin/data/format/hypfmt_e.html
%          for tremor, go to:
%           http://www-solid.eps.s.u-tokyo.ac.jp/~sloweq/?page=policy
%   
%       2. The time zone now used in both tremor and regular events are JST (UT+9),
%           Japan Standard Time = UTC + 9 h, which are consistent to each other, but need
%           conversion to UTC if individual earthquake analysises are needed
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/19
% Last modified date:   2020/02/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(2);

workpath = '/home/data2/chaosong/shikoku_kii';
tmrpath = strcat(workpath,'/etscat');
evtpath = strcat(workpath,'/jmacat02-18');
figpath = strcat(workpath,'/figs');

% region flag
regflag = 1;    % western shikoku
% regflag = 2;    % kii penninsula

if regflag == 1
    prefix = 'shikoku';
elseif regflag == 2
    prefix = 'kii';
end
disp(prefix);


% depth range that allows the regular earthquakes fall below the slab interface
depran = 10;

% recalculation flag
% recalflag = 1;
recalflag = 0;

% load the slab model
slab = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/ryu_slab1.0_clip.xyz');   % most similar


%% read regular event catalog
if recalflag    % if needs to recalculate the result
    % be careful about the time zone which are used
    evtall = Jmacatlogread(regflag,depran,slab);

else    % i.e. no need to recalculate, load the existing results
    % load the existing data
    data = load(strcat(evtpath,'/regevt',prefix,'dep',num2str(depran),'.mat'));
    evtall = data.evtall;
end
% ignore the ones with magnitude undetermined
evtall = evtall(evtall(:,11)~=1e6, :);

%% read tremor catalog from Obara NEID tremor
% from papers: 1. Obara, K., Tanaka, S., Maeda, T., & Matsuzawa, T. (2010)
% 2. Maeda, T., & Obara. K. (2009)
% IMPORTANT NOTICE
%   The time zone is JST (UT+9). I have confirmed that both catalogs are using the same
%	time zone by looking at the same detection 
if recalflag    % if needs to recalculate the result
    fpath = strcat(tmrpath,'/sloweq.20001231.6939.115617203.ObaraNEID.csv');
    obara = Sloweqread(fpath, regflag);
    
else
    data = load(strcat(tmrpath,'/tre',prefix,'.mat'));
    obara = data.treall;
end

%% make sure the regular event and tremor catalog have the same time coverage
% the ending time of regular events is later than tremor, so use the earlier one to be consistent
% if want to include ones with undetermined magnitude, ignore '& evtobj(:,11)~=1e6 '
evtall = evtall(evtall(:,1) >= obara(1,1) & evtall(:,1) <= obara(end,1), :);


%% plot tremor & regular events during ETS and obtain ones close in time & space
% set a empirical magnitude threshold from previous analysis 
minmag = 1.0;   % min magnitude
dttol = 2;      % max half time difference, in hr
dloctol = 0.05;  % max half location difference, in deg 

[stday, edday] = ETSestimate(regflag);
nets = size(stday,2);

evtets = [];
tmrets = [];
for i = 1: nets
% i=1;
    idbg = stday(i);      % iday of the beginning 
    ided = edday(i);
   
    f.fig=figure;
    widin = 6;  % maximum width allowed is 8.5 inches
    htin = 6;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/10 widin*res htin*res]);

    [f.ax,evtobj,tmrobj] = plttmrevt_spatiotemp(regflag,evtall,obara,idbg,ided,minmag,dttol,dloctol);
    evtets = [evtets; evtobj];
    tmrets = [tmrets; tmrobj];

end

%% find the events that are spatially close, but in tremor free days
% 'timexxx' for datetime type
% 'datexxx' for numerical date like yyyymmdd
% 'datestrxxx' for string type yyyymmdd
% 'idxxx' for the numerical order of one day relative to the first day
% 'ndxxx' for the number of days
time1 = datetime(obara(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
time2 = datetime(obara(end,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ndall = caldays(between(time1,time2,'days'))+1;

datetmr = unique(obara(:,1));
timetmr = datetime(datetmr,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
idtmr = caldays(between(time1,timetmr,'days'))+1;
idfree = setdiff(1:ndall,idtmr);
timefree = datetime(time1+caldays(idfree)-1,'Format','yyyy-MM-dd');
datefree = yyyymmdd(timefree');

minmag = 5;   % max magnitude difference
dloctol = 0.05;  % max half location difference, in deg
deptol = 2;     % max depth difference, in km
[evtobj] = evtouttmr(evtall, datefree, evtets, minmag, dloctol, deptol);



%% convert the time from local time zone to UTC
for i = 1: size(evtets,1)
    date = num2str(evtets(i,1));
    yr = date(1:4);
    mo = date(5:6);
    dy = date(7:8);
    if evtets(i,5) < 10
        hr = strcat('0',num2str(evtets(i,5)));
    else
        hr = num2str(evtets(i,5));
    end
    if evtets(i,6) < 10
        mi = strcat('0',num2str(evtets(i,6)));
    else
        mi = num2str(evtets(i,6));
    end
    sec = sprintf('%.3f',evtets(i,7));
    if evtets(i,7) < 10
        sec = strcat('0',sec);
    end
    timestr = strcat(yr,'-',mo,'-',dy,{' '},hr,':',mi,':',sec);

    d = datetime(timestr,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS','TimeZone','Asia/Tokyo');
    d.TimeZone = 'UTC'
%     dat(i) = d;
    
end




















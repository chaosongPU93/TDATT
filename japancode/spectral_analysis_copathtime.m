% function spectral_analysis_copathtime.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is to read the saved copath-time tremor and events; backgound
% tremor activity; other copath-events;
% Then it reads the corresponding sac file if i have it already, and convert
% it sac format if needed, remove instrument response, and finally calculate
% and plot their spectra;
% if the date is earlier than hinet online database, save these dates for
% distinct request methods;
% if the date is within hinet available query but i don't have it, then save
% the dates, then run python script 'down_hinet.py' to download them and then
% run the 'remove_response_downhinet.sh' to remove the station response for
% all data in the datelist. In order to keep tracking the record, a datelist
% for all dates that i don't have data previously and i use the new python
% script and bash shell. This is basically to not confuse myself, since the
% old data i have is the whole one-day 24 hr, but now to increase efficiency,
% i only download a short segment around the origin time, 10 min before, 20
% min after, 30min in total 
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/22
% Last modified date:   2020/07/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
clc;
% clear;
% close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = '/home/data2/chaosong/shikoku_kii';
tmrpath = strcat(workpath,'/etscat');
evtpath = strcat(workpath,'/jmacat02-18');
% figpath = strcat(workpath,'/figs');
figpath = strcat(workpath,'/target');
sacpath = '/home/data2/chaosong/japan_chao';
datapath = strcat(workpath,'/matsave');

% region flag
% regflag = 1;    % western shikoku
regflag = 2;    % kii penninsula

if regflag == 1
    prefix = 'shikoku';
elseif regflag == 2
    prefix = 'kii';
end
disp(prefix);


% depth range that allows the regular earthquakes fall below the slab interface
% depran = 15;
depran = inf;

% recalculation flag
% recalflag = 1;
recalflag = 0;


% load the slab model
slab = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/ryu_slab1.0_clip.xyz');   % most similar


%% read the regular events catalog
% 1. evtall = Jmacatlogread(regflag,depran,slab)
%   would read the events that are inside a rectangle defined by the 'regflag', and are within
%   a depth range below the slab interface defined by the 'depran' and 'slab'.
%   
% 2. [evtintmrall,evtouttmrall] = Jmacatlogread2(regflag,regbound,depran,slab)
%   even do more, you can divide the events into two exclusive parts by giving a boundary (now i am 
%   using the tremor active region as the boundary, but any closed boundary is allowed)

if recalflag    % if needs to recalculate the result
    % be careful about the time zone which are used
    [evtall,event] = Jmacatlogread(regflag,depran,slab);

else    % i.e. no need to recalculate, load the existing results
    % load the existing data
    data = load(strcat(evtpath,'/regevt',prefix,'dep',num2str(depran),'.mat'));
    evtall = data.evtall;
    
    data = load(strcat(evtpath,'/regevtstruct',prefix,'dep',num2str(depran),'.mat'));
    event = data.eventall;
end
% ignore the ones with magnitude undetermined
ind = find(evtall(:,11)~=1e6);
evtall = evtall(ind, :);
event = event(ind);

% read tremor catalog from Obara NEID tremor
if recalflag    % if needs to recalculate the result
    fpath = strcat(tmrpath,'/sloweq.20001231.6939.115617203.ObaraNEID.csv');
    obara = Sloweqread(fpath, regflag);
    
else
    data = load(strcat(tmrpath,'/tre',prefix,'.mat'));
    obara = data.treall;
end

% make sure the regular event and tremor catalog have the same time coverage
% the ending time of regular events is later than tremor, so use the earlier one to be consistent
% if want to include ones with undetermined magnitude,
ind = find(evtall(:,1) >= obara(1,1) & evtall(:,1) <= obara(end,1));
evtall = evtall(ind, :);
event = event(ind);

%%% check the event catalog more carefully, retain the natural EQ only
% event type, 1 for natural EQ, 5 for LFE
evttype = zeros(size(evtall,1),1);
for i = 1: size(evtall,1)
    evttype(i) = str2double(event(i).evttp);
end
evtall = evtall(evttype==1, :);     % 1 is natural EQ

% load all stations' information i have data
fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationlist'), 'r');
format = '%s %f %f %f \n';
catcell = textscan(fid, format);
fclose(fid);
stnm = char(catcell{1,1});    % station name, string
stlo = catcell{1,2};    % station longitude
stla = catcell{1,3};    % station latitude
stel = catcell{1,4};    % station elevation


%% load the target co-path-time tremor and event 
% velocity model used
velmod = 'jma2001';

% station chosen
sta = 'N.INMH';
% sta = 'N.HRKH';
% sta = 'N.HNZH';
% sta = 'N.TKWH';
% sta = 'N.KRTH';
% sta = 'N.KAWH';

for i = 1:size(stnm,1)
    if isequal(stnm(i,:), sta)
        stanum = i;
    end
end
stainfo{1,1} = stnm(stanum,:);
stainfo{1,2} = stlo(stanum);
stainfo{1,3} = stla(stanum);
stainfo{1,4} = stel(stanum);


% load saved results
rstmat = load(strcat(datapath,'/co-path-time_tmr_pierpt_',prefix,'_',sta,...
     '_',velmod,'.mat'));
evtintmr = rstmat.evtintmr;
tmrinevt = rstmat.tmrinevt;
tmrbkgd = rstmat.tmrbkgd;
evtbef = rstmat.evtbef;
evtaft = rstmat.evtaft;
evtbefgt = rstmat.evtbefgt;
evtaftgt = rstmat.evtaftgt;
evtbefle = rstmat.evtbefle;
evtaftle = rstmat.evtaftle;
evtbefgt10 = rstmat.evtbefgt10;
evtaftgt10 = rstmat.evtaftgt10;
evtbefle10 = rstmat.evtbefle10;
evtaftle10 = rstmat.evtaftle10;
clear rstmat


%% read the different types of event list
% read the complete date list for continuous data downloaded for kii
datedown = load(fullfile(sacpath,'trig_time_kii'));

% read the event list that i don't data before but i use down_hinet.py and
% remove_response_downhinet.sh to process it
rstmat = load(strcat(datapath,'/co-path-time_evtnodata',prefix,'_',sta,...
     '_',velmod,'.mat'));
evtnodata = rstmat.evtnodata;   % matrix for events i don't have data and need to run down_hinet.py to get them

evtbefhinet = [];   % matrix for events before hinet, i.e. earlier than 2004-04-01 

%% read data and obtain spectra
% set some parameters for data window
sps = 100;
lo = 0.5;
hi = 30;
npo = 2;
npa = 2;
% wlensec = 3.5;
% wlensecbf = 0.5;
wlensec = 5;
wlensecbf = 1;
pltflag = 1;
         
for i = 2: size(evtintmr,1)
    
    evtlabel = strcat(num2str(evtintmr(i,1)),num2str(evtintmr(i,5)), num2str(evtintmr(i,6)), ...
                      num2str(evtintmr(i,7)) );
                  
    if ismember(evtintmr(i,1),datedown)    % if event belongs to old data
        
        fprintf('%s, using existed old data.\n',evtlabel);
        
        %%% read one-day chunk of data
        % 120-s E, N, U components, p/P signal window, s/S signal window
        [eq1E,eq1N,eq1U,eq1t,eq1Eps,eq1Nps,eq1Ups,eq1tps,eq1Ess,eq1Nss,eq1Uss,eq1tss,...
            eq1Enoi,eq1Nnoi,eq1Unoi,eq1tnoi] = ...
            rdhinetevt(sacpath,stainfo,evtintmr(i,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag,...
            velmod);
        
    elseif ismember(evtintmr(i,:),evtnodata,'row')   % if event belongs to new data from down_hinet.py
        
        fprintf('%s, using newly downloaded data.\n',evtlabel);
        
        % read the event specific small segment of data
        % 120-s E, N, U components, p/P signal window, s/S signal window
        [eq1E,eq1N,eq1U,eq1t,eq1Eps,eq1Nps,eq1Ups,eq1tps,eq1Ess,eq1Nss,eq1Uss,eq1tss,...
            eq1Enoi,eq1Nnoi,eq1Unoi,eq1tnoi] = ...
            rdhinetevt_newdata(sacpath,stainfo,evtintmr(i,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,...
            pltflag,velmod);
        
    elseif ismember(evtintmr(i,1),datearch)
            
        fprintf('%s, using newly downloaded archive data.\n',evtlabel);
        
        % read the event specific small segment of data
        % 120-s E, N, U components, p/P signal window, s/S signal window
        [eq1E,eq1N,eq1U,eq1t,eq1Eps,eq1Nps,eq1Ups,eq1tps,eq1Ess,eq1Nss,eq1Uss,eq1tss,...
            eq1Enoi,eq1Nnoi,eq1Unoi,eq1tnoi] = ...
            rdhinetevt_newdata(sacpath,stainfo,evtintmr(i,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,...
            pltflag,velmod);
        
    end    

    % save figure if needed
    if pltflag
        print('-dpdf',strcat(figpath,'/Seismogram.',prefix,'_',sta,...
            '_',velmod,'_',evtlabel,'.pdf'),'-bestfit');
    end
    
    %%% obtain the PSD estimate and plot using pmtm and periodogram
    [pmeq1Ess,pmeq1Nss,pmeq1Uss,pmeq1Eps,pmeq1Nps,pmeq1Ups,pmft,f2,f3]=spectra_copathtime(eq1Eps,...
        eq1Nps,eq1Ups,eq1Ess,eq1Nss,eq1Uss,eq1Enoi,eq1Nnoi,eq1Unoi,sps); 
    
    % save figure if needed
    title(f2.ax(2),strcat(sta,{' periodgram '},evtlabel,{' mag='},num2str(evtintmr(i,11))));
    print(f2.fig,'-depsc2',strcat(figpath,'/PSD.pd.',prefix,'_',sta,...
        '_',velmod,'_',evtlabel,'.eps'));
    
    
    title(f3.ax(2),strcat(sta,{' pmtm '},evtlabel,{' mag='},num2str(evtintmr(i,11))));
    print(f3.fig,'-depsc2',strcat(figpath,'/PSD.pm.',prefix,'_',sta,...
        '_',velmod,'_',evtlabel,'.eps'));
    
    
end
         


           

















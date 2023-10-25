% function check_if_data_exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to read in the events that we are interested in, check
% if i already have data;
% if yes, convert the data format from win32 to sac, remove the station response
% and so on if necessary;
% if not, output those dates and events to run the corresponding python script
% to download the data specific to each event;
% Also, if the dates are earlier than hinet online database, then still output
% those dates, but downloading those data needs to send online query to hinet
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/24
% Last modified date:   2020/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

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


%% first, tackle the events that have co-path-time tremor
% read the complete date list for continuous data downloaded for kii
% this is the list that i already have data
datedown = load(fullfile(sacpath,'trig_time_kii'));

evtbefhinet = [];   % matrix for events before hinet, i.e. earlier than 2004-04-01
evtnodata = [];     % matrix for events i don't have data and need to run down_hinet.py to get them


% first, tackle the events that have co-path-time tremor
for i = 1: size(evtintmr,1)
    
    [evtbefhinet,evtnodata] = check_if_data_exist_func(evtintmr(i,:),sacpath,datedown,evtbefhinet,...
                                                       evtnodata);
                                                   
end

%% next, tackle the copath events before target event and closest tremor within 10 km is 5 hr before
for i = 2: size(evtbefgt10,1)
    tmpevt = evtbefgt10{i};
    for j = 1: size(tmpevt,1)
        
        [evtbefhinet,evtnodata] = check_if_data_exist_func(tmpevt(j,:),sacpath,datedown,evtbefhinet,...
                                                           evtnodata);
    end
end

%% next, tackle the copath events after target event and closest tremor within 10 km is 5 hr before
for i = 2: size(evtaftgt10,1)
    tmpevt = evtaftgt10{i};
    for j = 1: size(tmpevt,1)
        
        [evtbefhinet,evtnodata] = check_if_data_exist_func(tmpevt(j,:),sacpath,datedown,evtbefhinet,...
                                                           evtnodata);
    end
end



%% save data files
% output the event origin time for the events before 2004-04-01, precise to sec
fid = fopen(fullfile(sacpath,'datebefhinet'),'w+');
fprintf(fid,'%d %d %d %d %d %f \n',evtbefhinet(:,2:7)');
fclose(fid);

% output the event origin time for the events after 2004-04-01, but i don't have
% data, so i need to download them, precise to min
fid = fopen(fullfile(sacpath,'newdatelist'),'w+');
fprintf(fid,'%d %d %d %d %d \n',evtnodata(evtnodata(:,11)>=2,2:6)');
fclose(fid);

% output the event origin time for the events after 2004-04-01, but i don't have
% data, so i need to download them, precise to min
% append them to a combined date file as well
% fid = fopen(fullfile(sacpath,'datenewdownall'),'a+');
% fprintf(fid,'%d %d %d %d %d \n',evtnodata(:,2:6)');
% fclose(fid);

% save to mat file
save(strcat(datapath,'/co-path-time_evtbefhinet',prefix,'_',sta,...
    '_',velmod,'.mat'),'evtbefhinet');

save(strcat(datapath,'/co-path-time_evtnodata',prefix,'_',sta,...
    '_',velmod,'.mat'),'evtnodata');









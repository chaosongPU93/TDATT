% function spectral_analysis_copath.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is to read the saved copath-time tremor and events; backgound
% tremor activity; other copath-events and to obtain the spectra for copath
% events only.
% Like spectral_analysis_copathtime.m, It reads the corresponding sac file 
% if i have it already, calculate and plot their spectra;
%
% if the date is earlier than hinet online database, currently i haven't 
% found the way to do it, but should be not essentially different than others,
% since it just needs different method for downloading data
%
% if the date is within hinet available query i did't have, but i used
% python script 'down_hinet.py' to download them and then read them in a 
% slightly different way, thne calculate and plot their spectra;
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/25
% Last modified date:   2020/07/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
% clear;
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


%% read the corresponding event list if i have it downloaded
% read the complete date list for continuous data downloaded for kii
datedown = load(fullfile(sacpath,'trig_time_kii'));

% read the event list that i don't data before but i use down_hinet.py and
% remove_response_downhinet.sh to process it
rstmat = load(strcat(datapath,'/co-path-time_evtnodata',prefix,'_',sta,...
     '_',velmod,'.mat'));
evtnodata = rstmat.evtnodata;   % matrix for events i don't have data and need to run down_hinet.py to get them

% evtbefhinet = [];   % matrix for events before hinet, i.e. earlier than 2004-04-01 

% read the date list of archive data from hinet before 20040401
datearch = load(fullfile(sacpath,'archlist'));

%% 
% set some parameters for data window
sps = 100;
lo = 0.5;
hi = 30;
npo = 2;
npa = 2;
wlensec = 5;
wlensecbf = 1;
pltflag = 0;

%% first tackle the co-path events that occurred before the target event
eq1 = cell(size(evtbefgt10,1),16);
for i = 2: size(evtbefgt10,1)
    tmpevt = evtbefgt10{i};
    
%     for j = 1: size(tmpevt,1)
        j=1;
        evtlabel = strcat(num2str(tmpevt(j,1)),num2str(tmpevt(j,5)), num2str(tmpevt(j,6)), ...
                          num2str(tmpevt(j,7)) );
        disp(tmpevt(j,1))
        if ismember(tmpevt(j,1),datedown)    % if event belongs to old data
            
            fprintf('%s, using existed old data.\n',evtlabel);
            
            %%% read one-day chunk of data
            % 120-s E, N, U components, p/P signal window, s/S signal window
            [eq1E,eq1N,eq1U,eq1t,eq1Eps,eq1Nps,eq1Ups,eq1tps,eq1Ess,eq1Nss,eq1Uss,eq1tss,...
                eq1Enoi,eq1Nnoi,eq1Unoi,eq1tnoi] = ...
                rdhinetevt(sacpath,stainfo,tmpevt(j,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag,...
                velmod);
            
        elseif ismember(tmpevt(j,:),evtnodata,'row')   % if event belongs to new data from down_hinet.py
                    
            fprintf('%s, using newly downloaded data.\n',evtlabel);
        
            % read the event specific small segment of data
            % 120-s E, N, U components, p/P signal window, s/S signal window
            [eq1E,eq1N,eq1U,eq1t,eq1Eps,eq1Nps,eq1Ups,eq1tps,eq1Ess,eq1Nss,eq1Uss,eq1tss,...
                eq1Enoi,eq1Nnoi,eq1Unoi,eq1tnoi] = ...
                rdhinetevt_newdata(sacpath,stainfo,tmpevt(j,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,...
                pltflag,velmod);
        
        end

        % combine all results into one matrix    
        eq1EA(:,j) = eq1E;
        eq1NA(:,j) = eq1N;
        eq1UA(:,j) = eq1U;
        eq1tA(:,j) = eq1t;
        eq1EpsA(:,j) = eq1Eps;
        eq1NpsA(:,j) = eq1Nps;
        eq1UpsA(:,j) = eq1Ups;
        eq1tpsA(:,j) = eq1tps;
        eq1EssA(:,j) = eq1Ess;
        eq1NssA(:,j) = eq1Nss;
        eq1UssA(:,j) = eq1Uss;
        eq1tssA(:,j) = eq1tss;
        eq1EnoiA(:,j) = eq1Enoi;
        eq1NnoiA(:,j) = eq1Nnoi;
        eq1UnoiA(:,j) = eq1Unoi;
        eq1tnoiA(:,j) = eq1tnoi;

        % save figure if needed
%         if pltflag
%             print('-dpdf',strcat(figpath,'/Seismogram.',prefix,'_',sta,...
%                   '_',velmod,'_',evtlabel,'.pdf'),'-bestfit');
%         end
    
%     end    
    
    %%% obtain the PSD estimate and plot using pmtm and periodogram
    [pmeq1EssA,pmeq1NssA,pmeq1UssA,pmeq1EpsA,pmeq1NpsA,pmeq1UpsA,pmftA,f1,f2]=spectra_copath(eq1EpsA,...
            eq1NpsA,eq1UpsA,eq1EssA,eq1NssA,eq1UssA,eq1EnoiA,eq1NnoiA,eq1UnoiA,sps,tmpevt);
    
    % save figure if needed
    supertit(f1.ax,strcat(sta,{' periodgram '},evtlabel));
    print(f1.fig,'-depsc2',strcat(figpath,'/PSD.pd.',prefix,'_',sta,...
          '_',velmod,'_',evtlabel,'.eps'));
    
    
    supertit(f2.ax,strcat(sta,{' pmtm '},evtlabel));
    print(f2.fig,'-depsc2',strcat(figpath,'/PSD.pm.',prefix,'_',sta,...
              '_',velmod,'_',evtlabel,'.eps'));
    
    
end
         

%% Then tackle the co-path events that occurred after the target event
eq2 = cell(size(evtaftgt10,1),16);
for i = 2: size(evtaftgt10,1)
    tmpevt = evtaftgt10{i};
    
    tmpevt = tmpevt(tmpevt(:,11)>=2, :);
    
    eq2EA = [];
    eq2NA = [];
    eq2UA = [];
    eq2tA = [];
    eq2EpsA = [];
    eq2NpsA = [];
    eq2UpsA = [];
    eq2tpsA = [];
    eq2EssA = [];
    eq2NssA = [];
    eq2UssA = [];
    eq2tssA = [];
    eq2EnoiA = [];
    eq2NnoiA = [];
    eq2UnoiA = [];
    eq2tnoiA = [];
    noemptyi = [];   % index of normal event
    
    for j = 1: size(tmpevt,1)
%         j=2;
        evtlabel = strcat(num2str(tmpevt(j,1)),num2str(tmpevt(j,5)), num2str(tmpevt(j,6)), ...
                          num2str(tmpevt(j,7)) );
        if ismember(tmpevt(j,1),datedown)    % if event belongs to old data
            
            fprintf('%s, using existed old data.\n',evtlabel);
            
            %%% read one-day chunk of data
            % 120-s E, N, U components, p/P signal window, s/S signal window
            [eq2E,eq2N,eq2U,eq2t,eq2Eps,eq2Nps,eq2Ups,eq2tps,eq2Ess,eq2Nss,eq2Uss,eq2tss,...
                eq2Enoi,eq2Nnoi,eq2Unoi,eq2tnoi] = ...
                rdhinetevt(sacpath,stainfo,tmpevt(j,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag,...
                velmod);
            
        elseif ismember(tmpevt(j,:),evtnodata,'row')   % if event belongs to new data from down_hinet.py
                    
            fprintf('%s, using newly downloaded data.\n',evtlabel);
        
            % read the event specific small segment of data
            % 120-s E, N, U components, p/P signal window, s/S signal window
            [eq2E,eq2N,eq2U,eq2t,eq2Eps,eq2Nps,eq2Ups,eq2tps,eq2Ess,eq2Nss,eq2Uss,eq2tss,...
                eq2Enoi,eq2Nnoi,eq2Unoi,eq2tnoi] = ...
                rdhinetevt_newdata(sacpath,stainfo,tmpevt(j,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,...
                pltflag,velmod);
        elseif ismember(tmpevt(j,1),datearch)
            
            fprintf('%s, using newly downloaded archive data.\n',evtlabel);
            
            % read the event specific small segment of data
            % 120-s E, N, U components, p/P signal window, s/S signal window
            [eq2E,eq2N,eq2U,eq2t,eq2Eps,eq2Nps,eq2Ups,eq2tps,eq2Ess,eq2Nss,eq2Uss,eq2tss,...
                eq2Enoi,eq2Nnoi,eq2Unoi,eq2tnoi] = ...
                rdhinetevt_newdata(sacpath,stainfo,tmpevt(j,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,...
                pltflag,velmod);
            
%             fprintf('%s, NO DATA.\n',evtlabel);
%             eq2E = [];
%             eq2N = [];
%             eq2U = [];
%             eq2t = [];
%             eq2Eps = [];
%             eq2Nps = [];
%             eq2Ups = [];
%             eq2tps = [];
%             eq2Ess = [];
%             eq2Nss = [];
%             eq2Uss = [];
%             eq2tss = [];
%             eq2Enoi = [];
%             eq2Nnoi = [];
%             eq2Unoi = [];
%             eq2tnoi = [];
        
        end

        % combine all results into one matrix
%         if ~isempty(eq1E) && ~isempty(eq1Eps) && ~isempty(eq1Ess) && ~isempty(eq1Enoi)
%         eq1EA(:,j) = eq1E;
%         eq1NA(:,j) = eq1N;
%         eq1UA(:,j) = eq1U;
%         eq1tA(:,j) = eq1t;
%         eq1EpsA(:,j) = eq1Eps;
%         eq1NpsA(:,j) = eq1Nps;
%         eq1UpsA(:,j) = eq1Ups;
%         eq1tpsA(:,j) = eq1tps;
%         eq1EssA(:,j) = eq1Ess;
%         eq1NssA(:,j) = eq1Nss;
%         eq1UssA(:,j) = eq1Uss;
%         eq1tssA(:,j) = eq1tss;
%         eq1EnoiA(:,j) = eq1Enoi;
%         eq1NnoiA(:,j) = eq1Nnoi;
%         eq1UnoiA(:,j) = eq1Unoi;
%         eq1tnoiA(:,j) = eq1tnoi;
%          end
        eq2EA = [eq2EA eq2E];
        eq2NA = [eq2NA eq2N];
        eq2UA = [eq2UA eq2U];
        eq2tA = [eq2tA eq2t];
        eq2EpsA = [eq2EpsA eq2Eps];
        eq2NpsA = [eq2NpsA eq2Nps];
        eq2UpsA = [eq2UpsA eq2Ups];
        eq2tpsA = [eq2tpsA eq2tps];
        eq2EssA = [eq2EssA eq2Ess];
        eq2NssA = [eq2NssA eq2Nss];
        eq2UssA = [eq2UssA eq2Uss];
        eq2tssA = [eq2tssA eq2tss];
        eq2EnoiA = [eq2EnoiA eq2Enoi];
        eq2NnoiA = [eq2NnoiA eq2Nnoi];
        eq2UnoiA = [eq2UnoiA eq2Unoi];
        eq2tnoiA = [eq2tnoiA eq2tnoi];

        if length(eq2Ess) == wlensec*sps
            noemptyi = [noemptyi; j]; 
        end
        % save figure if needed
%         if pltflag
%             print('-dpdf',strcat(figpath,'/Seismogram.',prefix,'_',sta,...
%                   '_',velmod,'_',evtlabel,'.pdf'),'-bestfit');
%         end
    
    end    
    
    %%% obtain the PSD estimate and plot using pmtm and periodogram
    [pmeq2EssA,pmeq2NssA,pmeq2UssA,pmeq2EpsA,pmeq2NpsA,pmeq2UpsA,pmft,f2,f3,f4]=spectra_copath(eq2EpsA,...
            eq2NpsA,eq2UpsA,eq2EssA,eq2NssA,eq2UssA,eq2EnoiA,eq2NnoiA,eq2UnoiA,sps,tmpevt(noemptyi,:));
    
    % save figure if needed
    evtlblmas = strcat(num2str(evtintmr(i,1)),num2str(evtintmr(i,5)), num2str(evtintmr(i,6)), ...
                       num2str(evtintmr(i,7)) );
    
%     title(f2.ax(2),strcat(sta,{' pmtm--'},{'evts after '},evtlblmas,{' p1'}));
%     print(f2.fig,'-depsc2',strcat(figpath,'/PSD.pm.p1',prefix,'_',sta,...
%           '_',velmod,'_evtsaft_',evtlblmas,'.eps'));
    

    title(f3.ax(2),strcat(sta,{' pmtm--'},{'evts after '},evtlblmas,{' p2'}));
    print(f3.fig,'-depsc2',strcat(figpath,'/PSD.pm.p2',prefix,'_',sta,...
          '_',velmod,'_evtsaft_',evtlblmas,'.eps'));
      
%     title(f4.ax(2),strcat(sta,{' pmtm--'},{'evts after '},evtlblmas,{' p3'}));
%     print(f4.fig,'-depsc2',strcat(figpath,'/PSD.pm.p3',prefix,'_',sta,...
%           '_',velmod,'_evtsaft_',evtlblmas,'.eps'));  
    
    
end
           

















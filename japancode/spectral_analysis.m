% function spectral_analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to read the data from Japan Hinet (after or before instrument
% response removal), choose segments of records according to the origin time of 
% regular earthquakes or tremor detections, and analyze and compare their spectra.
% Hopefully we could find some statistics on difference or similarity behind these
% spectra
%
% 
%   NOTICE:
%       1. The time zone now used in both tremor, regular events catalog and the Hinet
%           data are JST (UT+9), Japan Standard Time = UTC + 9 h, which are consistent 
%           to each other. Unless specified, conversion to UTC frame is unnecessary
%
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/03/03
% Last modified date:   2020/03/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = '/home/data2/chaosong/shikoku_kii';
tmrpath = strcat(workpath,'/etscat');
evtpath = strcat(workpath,'/jmacat02-18');
figpath = strcat(workpath,'/figs');
sacpath = '/home/data2/chaosong/japan_chao';

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
depran = 10;

% recalculation flag
% recalflag = 1;
recalflag = 0;

% load the slab model
slab = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/ryu_slab1.0_clip.xyz');   % most similar

% load chosen station names only, not all of them
fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationlistsel'), 'r');
format = '%s %f %f %f \n';
stainfo = textscan(fid, format);
fclose(fid);
stnmsel = char(stainfo{1,1});    % station name, string
stlosel = stainfo{1,2};    % station longitude
stlasel = stainfo{1,3};    % station latitude
stelsel = stainfo{1,4};    % station elevation

% read tremor catalog from Obara NEID tremor
if recalflag    % if needs to recalculate the result
    fpath = strcat(tmrpath,'/sloweq.20001231.6939.115617203.ObaraNEID.csv');
    obara = Sloweqread(fpath, regflag);
    
else
    data = load(strcat(tmrpath,'/tre',prefix,'.mat'));
    obara = data.treall;
end

% get regular events vertically within the depran below, horizontally within/out the tremor site
if recalflag    % if needs to recalculate the result
    % be careful about the time zone which are used
    % evtinall is the events inside the tremor region, whereas evtoutall is those outside  
    [evtinall,evtoutall] = Jmacatlogread2(regflag,tmrbd,depran,slab);

else    % i.e. no need to recalculate, load the existing results
    % load the existing data
    data = load(strcat(evtpath,'/regevt_in_tmrbound',prefix,'dep',num2str(depran),'.mat'));
    evtinall = data.evtintmrall;
    
    data = load(strcat(evtpath,'/regevt_out_tmrbound',prefix,'dep',num2str(depran),'.mat'));
    evtoutall = data.evtouttmrall;

end
% ignore the ones with magnitude undetermined
evtinall = evtinall(evtinall(:,11)~=1e6, :);
evtoutall = evtoutall(evtoutall(:,11)~=1e6, :);

% make sure the regular event and tremor catalog have the same time coverage
% the ending time of regular events is later than tremor, so use the earlier one to be consistent
% if want to include ones with undetermined magnitude, 
evtinall = evtinall(evtinall(:,1) >= obara(1,1) & evtinall(:,1) <= obara(end,1), :);
evtoutall = evtoutall(evtoutall(:,1) >= obara(1,1) & evtoutall(:,1) <= obara(end,1), :);

% evtall is all regular events in the rectangle region 
evtall = [evtinall; evtoutall];
evtall = sortrows(evtall,1);


% load objective events for different purposes
data = load(strcat(evtpath,'/objevts_forspectra.',prefix,'dep',num2str(depran),'.mat'));
evtir = data.evtir;     % events that happen close to a tremor (+-0.01 deg) within +-2 hr
% evtir = sortrows(evtir,1);
evts = data.evts;       % events that are spatially close, but in tremor free days
% evts = sortrows(evts,1);
evtti = data.evtti;     % events that are timely close, but outside the in radius regions
% evtti = sortrows(evtti,1);
evtto = data.evtto;     % events that are timely close, but outside the entire tremor regions
% evtto = sortrows(evtto,1);

minmag = 1;
evtirm1 = evtir(evtir(:,11)>=minmag, :);  % events that happen close to a tremor (+-0.01 deg) within +-2 hr
% evtirm1 = sortrows(evtirm1,1);
evtsm1 = evts(evts(:,11)>=minmag, :);   % events that are spatially close, but in tremor free days
% evtsm1 = sortrows(evtsm1,1);
evttim1 = evtti(evtti(:,11)>=minmag, :); % events that are timely close, but outside the in radius regions
% evttim1 = sortrows(evttim1,1);
evttom1 = evtto(evtto(:,11)>=minmag, :); % events that are timely close, but outside the entire tremor regions
% evttom1 = sortrows(evttom1,1);


dtevtirm1 = unique(evtirm1(:,1));
dtevttom1 = unique(evttom1(:,1));
dtevtsm1 = unique(evtsm1(:,1));
dtevttim1 = unique(evttim1(:,1));

% read the complete date list for continuous data downloaded for kii
datedown = load(fullfile(sacpath,'trig_time_kii'));

evtirava = [];
evtsava = [];
% evttiava = [];
evttoava = [];
for i = 1:size(evtir,1)
    tmp1 = evts(evts(:,end)==evtir(i,end), :);
%     tmp2 = evtti(evtti(:,end)==evtir(i,end), :); && ~isempty(tmp2)
    tmp3 = evtto(evtto(:,end)==evtir(i,end), :);
    if ~isempty(tmp1) &&  ~isempty(tmp3)
        evtirava = [evtirava; evtir(i,:)];
        evtsava = [evtsava; tmp1];
%         evttiava = [evttiava; tmp2];
        evttoava = [evttoava; tmp3];
    end
end

% evtirtest = evtirava(evtirava(:,11)==max(evtirava(:,11)), :);
% 
% evtstest = evtsava(abs(evtsava(:,11)-evtirtest(:,11)) <= 0.1 & ...
%                    abs(evtsava(:,10)-evtirtest(:,10)) <= 2 &...
%                    evtsava(:,end)==evtirtest(:,end), :);
%                
% evttotest = evttoava(abs(evttoava(:,11)-evtirtest(:,11)) <= 0.1 & ...
%                    abs(evttoava(:,10)-evtirtest(:,10)) <= 5 &...
%                    evttoava(:,end)==evtirtest(:,end), :);
               
               
evtirtest = evtirava(end,:);

evtstest = evtsava(abs(evtsava(:,11)-evtirtest(:,11)) <= 0.2 & ...
                   abs(evtsava(:,10)-evtirtest(:,10)) <= 2 &...
                   evtsava(:,end)==evtirtest(:,end), :);
               
evttotest = evttoava(abs(evttoava(:,11)-evtirtest(:,11)) <= 1 & ...
                   abs(evttoava(:,10)-evtirtest(:,10)) <= 5 &...
                   evttoava(:,end)==evtirtest(:,end), :);
               



%% evtir or evtirm1, regular EQ1 in tremor days: tremor + earthquake + noise
% events that happen close to a tremor (+-0.01 deg) within +-2 hr

% dateava = intersect(evtirtest(:,1),datedown);

sps = 50;
lo = 0.05;
hi = 20;
npo = 2;
npa = 2;
wlensec = 2.5;
wlensecbf = 0.5;
pltflag = 1;

% earthquake type 1, 120-s E, N, U components, p/S signal window, s/S signal window
[eq1E,eq1N,eq1U,eq1t,eq1Eps,eq1Nps,eq1Ups,eq1tps,eq1Ess,eq1Nss,eq1Uss,eq1tss,...
        eq1Enoi,eq1Nnoi,eq1Unoi,eq1tnoi] = ...
    rdhinetevt(sacpath,stainfo,evtirtest,sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag);


%% tremor of the same day: tremor - earthquake + noise, sequential average

sps = 50;
lo = 0.05;
hi = 20;
npo = 2;
npa = 2;
periodsec = 400;
pltflag = 1;

[tmrE,tmrN,tmrU,tmrt] = ...
    rdhinettmr(sacpath,evtirtest(1),stainfo,evtall,obara,sps,lo,hi,npo,npa,periodsec,pltflag);



%% evts or evtsm1, regular EQ 2 in tremor free days: -tremor + earthquake + noise
% events that are spatially close, but in tremor free days

% evttmp = evtsm1(evtsm1(:,end)==evtirm1(1,end), :);
% 
% dateava = intersect(evttmp(:,1),datedown);
% if ~isempty(dateava)
%     evtava = [];
%     for i = 1: length(dateava)
%         tmp = evttmp(evttmp(:,1)==dateava(i), :);
%         evtava = [evtava; tmp];
%     end
% 
%     datetest = evtava(1,1);
%     
% end

sps = 50;
lo = 0.05;
hi = 20;
npo = 2;
npa = 2;
wlensec = 2.5;
wlensecbf = 0.5;
pltflag = 1;

% earthquake type 2, 120-s E, N, U components, p/S signal window, s/S signal window
[eq2E,eq2N,eq2U,eq2t,eq2Eps,eq2Nps,eq2Ups,eq2tps,eq2Ess,eq2Nss,eq2Uss,eq2tss,...
        eq2Enoi,eq2Nnoi,eq2Unoi,eq2tnoi] = ...
    rdhinetevt(sacpath,stainfo,evtstest(2,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag);



% %% evtti or evttim1, regular EQ 3 in tremor free days: -tremor + earthquake + noise
% % events that are timely close, but outside the in radius regions
% 
% evttmp = evtti(evtti(:,end)==evtirm1(1,end), :);
% evttmp = evttmp(evttmp(:,11)==max(evttmp(:,11)), :);
% 
% dateava = intersect(evttmp(:,1),datedown);
% if ~isempty(dateava)
%     evtava = [];
%     for i = 1: length(dateava)
%         tmp = evttmp(evttmp(:,1)==dateava(i), :);
%         evtava = [evtava; tmp];
%     end
% 
%     datetest = evtava(1,1);
%     
% end
% 
% 


%% evtto or evttom1, regular EQ 4 in tremor free days: -tremor + earthquake + noise
% events that are timely close, but outside the entire tremor regions

sps = 50;
lo = 0.05;
hi = 20;
npo = 2;
npa = 2;
wlensec = 2.5;
wlensecbf = 0.5;
pltflag = 1;

% earthquake type 2, 120-s E, N, U components, p/S signal window, s/S signal window
[eq4E,eq4N,eq4U,eq4t,eq4Eps,eq4Nps,eq4Ups,eq4tps,eq4Ess,eq4Nss,eq4Uss,eq4tss,...
          eq4Enoi,eq4Nnoi,eq4Unoi,eq4tnoi] = ...
    rdhinetevt(sacpath,stainfo,evttotest(1,:),sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag);


%% ambient noise in tremor free days:   -tremor - earthquake + noise
d1 = datetime(obara(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
d2 = datetime(obara(end,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
nday = caldays(between(d1,d2,'days'))+1; 
iallday = 1:nday;
iallday = iallday';

% day sequence of all tremor detections
datetmr = unique(obara(:,1));
dtmr = datetime(datetmr,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
itmrday = caldays(between(d1,dtmr,'days'))+1; 

% day sequence of all events
dateevt = unique(evtall(:,1));
devt = datetime(dateevt,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ievtday = caldays(between(d1,devt,'days'))+1; 

ievtmday = union(itmrday,ievtday);  % days have either tremor or events
ifreeday = setdiff(iallday, ievtmday);  % days have neither tremor nor events

dfree = datetime(d1+caldays(ifreeday)-1,'Format','yyyy-MM-dd');
datefree = yyyymmdd(dfree);     % dates in numbers yyyymmdd


% find their intersections, if empty, then new data needs to be downloaded 
dateava = intersect(datefree,datedown);

if isempty(dateava)
    disp('Have no data for tremor&event days, need to download new data')    
end

    
datetest = dateava(2);

sps = 50;
lo = 0.05;
hi = 20;
npo = 2;
npa = 2;
pltflag = 1;

[noiEa,noiNa,noiUa,noita] = rdhinetnoise(sacpath,stainfo,datetest,sps,lo,hi,npo,npa,pltflag);
sttime = 8e4;
dt = 1/sps;
stind = find(abs(noita(:,1)-sttime)<=dt,1,'first');
periodsec = 400;
noiE = noiEa(stind: stind+periodsec*sps-1, :);
noiN = noiNa(stind: stind+periodsec*sps-1, :);
noiU = noiUa(stind: stind+periodsec*sps-1, :);
noit = noita(stind: stind+periodsec*sps-1, :);




%% obtain the spectra of different earthquakes, tremor and noise via periodgram
% For PSD, use pmtm or periodgram (no averaging) for regular earthquakes; use pchave for tremor and
% noise;
% For amplitude spectrum, use fft, see https://www.mathworks.com/help/matlab/ref/fft.html

%%%% PSD estimate with periodgram
%%% EQ 1
% periodogram
nfft = pow2(nextpow2(size(eq1Ess,1))-1);
window = hann(size(eq1Ess,1));
Fs = sps;
[pdeq1Ess,pdft] = periodogram(eq1Ess,window,nfft,Fs);  % psd is done to each col, i.e., each station 
[pdeq1Nss,~] = periodogram(eq1Nss,window,nfft,Fs);
[pdeq1Uss,~] = periodogram(eq1Uss,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq1Ess,pmft] = pmtm(eq1Ess,nw,nfft,Fs);
[pmeq1Nss,~] = pmtm(eq1Nss,nw,nfft,Fs);
[pmeq1Uss,~] = pmtm(eq1Uss,nw,nfft,Fs);


%%% EQ 2
% periodogram
nfft = pow2(nextpow2(size(eq2Ess,1))-1);
window = hann(size(eq2Ess,1));
Fs = sps;
[pdeq2Ess,~] = periodogram(eq2Ess,window,nfft,Fs);
[pdeq2Nss,~] = periodogram(eq2Nss,window,nfft,Fs);
[pdeq2Uss,~] = periodogram(eq2Uss,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq2Ess,~] = pmtm(eq2Ess,nw,nfft,Fs);
[pmeq2Nss,~] = pmtm(eq2Nss,nw,nfft,Fs);
[pmeq2Uss,~] = pmtm(eq2Uss,nw,nfft,Fs);


%%% EQ 4
% periodogram
nfft = pow2(nextpow2(size(eq4Ess,1))-1);
window = hann(size(eq4Ess,1));
Fs = sps;
[pdeq4Ess,~] = periodogram(eq4Ess,window,nfft,Fs);
[pdeq4Nss,~] = periodogram(eq4Nss,window,nfft,Fs);
[pdeq4Uss,~] = periodogram(eq4Uss,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq4Ess,~] = pmtm(eq4Ess,nw,nfft,Fs);
[pmeq4Nss,~] = pmtm(eq4Nss,nw,nfft,Fs);
[pmeq4Uss,~] = pmtm(eq4Uss,nw,nfft,Fs);


%%% EQ 1, noise around origin
% periodogram
nfft = pow2(nextpow2(size(eq1Enoi,1))-1);
window = hann(size(eq1Enoi,1));
Fs = sps;
[pdeq1Enoi,~] = periodogram(eq1Enoi,window,nfft,Fs);
[pdeq1Nnoi,~] = periodogram(eq1Nnoi,window,nfft,Fs);
[pdeq1Unoi,~] = periodogram(eq1Unoi,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq1Enoi,~] = pmtm(eq1Enoi,nw,nfft,Fs);
[pmeq1Nnoi,~] = pmtm(eq1Nnoi,nw,nfft,Fs);
[pmeq1Unoi,~] = pmtm(eq1Unoi,nw,nfft,Fs);


%%% EQ 2, noise around origin
% periodogram
nfft = pow2(nextpow2(size(eq2Enoi,1))-1);
window = hann(size(eq2Enoi,1));
Fs = sps;
[pdeq2Enoi,~] = periodogram(eq2Enoi,window,nfft,Fs);
[pdeq2Nnoi,~] = periodogram(eq2Nnoi,window,nfft,Fs);
[pdeq2Unoi,~] = periodogram(eq2Unoi,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq2Enoi,~] = pmtm(eq2Enoi,nw,nfft,Fs);
[pmeq2Nnoi,~] = pmtm(eq2Nnoi,nw,nfft,Fs);
[pmeq2Unoi,~] = pmtm(eq2Unoi,nw,nfft,Fs);


%%% EQ 4, noise around origin
% periodogram
nfft = pow2(nextpow2(size(eq4Enoi,1))-1);
window = hann(size(eq4Enoi,1));
Fs = sps;
[pdeq4Enoi,~] = periodogram(eq4Enoi,window,nfft,Fs);
[pdeq4Nnoi,~] = periodogram(eq4Nnoi,window,nfft,Fs);
[pdeq4Unoi,~] = periodogram(eq4Unoi,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq4Enoi,~] = pmtm(eq4Enoi,nw,nfft,Fs);
[pmeq4Nnoi,~] = pmtm(eq4Nnoi,nw,nfft,Fs);
[pmeq4Unoi,~] = pmtm(eq4Unoi,nw,nfft,Fs);


%%% tremor, mean of the consecutive 2.5-s window over a 400-s time period
% periodogram
wlensec = 2.5;
wlen = wlensec*sps;
nwin = round(size(tmrE,1)/wlen);

nfft = pow2(nextpow2(wlen)-1);
window = hann(wlen);
Fs = sps;
pdtmrE = zeros(round(nfft/2+1),size(tmrE,2));
pdtmrN = zeros(round(nfft/2+1),size(tmrE,2));
pdtmrU = zeros(round(nfft/2+1),size(tmrE,2));
for i = 1: nwin
    tmpe = tmrE((i-1)*wlen+1: i*wlen, :);
    tmpn = tmrN((i-1)*wlen+1: i*wlen, :);
    tmpu = tmrU((i-1)*wlen+1: i*wlen, :);
    
    [ptmpe,~] = periodogram(tmpe,window,nfft,Fs);
    [ptmpn,~] = periodogram(tmpn,window,nfft,Fs);
    [ptmpu,~] = periodogram(tmpu,window,nfft,Fs);
    
    pdtmrE = pdtmrE+ptmpe;
    pdtmrN = pdtmrN+ptmpe;
    pdtmrU = pdtmrU+ptmpe;
    
end
pdtmrE = pdtmrE./nwin;
pdtmrN = pdtmrN./nwin;
pdtmrU = pdtmrU./nwin;

% pchave
olap = 50;
pctmrE = zeros(round(nfft/2+1),size(tmrE,2));
pctmrN = zeros(round(nfft/2+1),size(tmrE,2));
pctmrU = zeros(round(nfft/2+1),size(tmrE,2));
for i = 1: size(tmrE,2)
    [pctmrE(:,i),pcft]=pchave(tmrE(:,i),wlen,olap,nfft,Fs,'MAD','dpss');
    [pctmrN(:,i),~]=pchave(tmrN(:,i),wlen,olap,nfft,Fs,'MAD','dpss');
    [pctmrU(:,i),~]=pchave(tmrU(:,i),wlen,olap,nfft,Fs,'MAD','dpss');
end


%%% ambient noise, mean of the consecutive 2.5-s window over a 400-s time period
% periodogram
wlensec = 2.5;
wlen = wlensec*sps;
nwin = round(size(noiE,1)/wlen);

nfft = pow2(nextpow2(wlen)-1);
window = hann(wlen);
Fs = sps;
pdnoiE = zeros(round(nfft/2+1),size(noiE,2));
pdnoiN = zeros(round(nfft/2+1),size(noiE,2));
pdnoiU = zeros(round(nfft/2+1),size(noiE,2));
for i = 1: nwin
    tmpe = noiE((i-1)*wlen+1: i*wlen, :);
    tmpn = noiN((i-1)*wlen+1: i*wlen, :);
    tmpu = noiU((i-1)*wlen+1: i*wlen, :);
    
    [ptmpe,~] = periodogram(tmpe,window,nfft,Fs);
    [ptmpn,~] = periodogram(tmpn,window,nfft,Fs);
    [ptmpu,~] = periodogram(tmpu,window,nfft,Fs);
    
    pdnoiE = pdnoiE+ptmpe;
    pdnoiN = pdnoiN+ptmpe;
    pdnoiU = pdnoiU+ptmpe;
    
end
pdnoiE = pdnoiE./nwin;
pdnoiN = pdnoiN./nwin;
pdnoiU = pdnoiU./nwin;

% pchave
olap = 50;
pcnoiE = zeros(round(nfft/2+1),size(noiE,2));
pcnoiN = zeros(round(nfft/2+1),size(noiE,2));
pcnoiU = zeros(round(nfft/2+1),size(noiE,2));
for i = 1: size(tmrE,2)
    [pcnoiE(:,i),pcft]=pchave(noiE(:,i),wlen,olap,nfft,Fs,'MAD','dpss');
    [pcnoiN(:,i),~]=pchave(noiN(:,i),wlen,olap,nfft,Fs,'MAD','dpss');
    [pcnoiU(:,i),~]=pchave(noiU(:,i),wlen,olap,nfft,Fs,'MAD','dpss');
end


%% plot periodgram PSD 

ifig = 0;
for i = 1: size(stnmsel,1)     % each station
    
    if mod(i,3)==1
       
        f.fig=figure;
        f.fig.Renderer='Painters';
        widin = 8.3;  % maximum width allowed is 8.5 inches
        htin = 9;   % maximum height allowed is 11 inches
        set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

        ifig = ifig+1;
    end
    
    nrow = 3;     % number of stations
    ncol = 3;   % number of components
    
    %%% N component
    
    j=1;
    isub = (i-(ifig-1)*nrow-1)*ncol+j;
    subplot(nrow,ncol,isub)
    box on
    % S signal
    semilogx(pdft,10*log10(pdeq1Nss(:,i)),'-','color','k','linew',1); hold on
    semilogx(pdft,10*log10(pdeq2Nss(:,i)),'-','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pdft,10*log10(pdeq4Nss(:,i)),'-','color',[0.7 0.7 0.7],'linew',1);
    % S noise before origin
    semilogx(pdft,10*log10(pdeq1Nnoi(:,i)),'-.','color','k','linew',1);
    semilogx(pdft,10*log10(pdeq2Nnoi(:,i)),'-.','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pdft,10*log10(pdeq4Nnoi(:,i)),'-.','color',[0.7 0.7 0.7],'linew',1);
    % tremor
    semilogx(pdft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
    % ambient noise
    semilogx(pdft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);
    
    xticks([1e0, 1e1]);
    xlim([0.5, Fs/2]);
%     ylim(f.ax(isub),[-150,-20]);
%     xlabel(f.ax(isub),'Frequency (Hz)');
%     ylabel(f.ax(isub),'PSD (dB/Hz)');
    title(strcat(stnmsel(i,:),{'.N periodgram'}));
    lgd = legend({'evtir','evts','evtto','noiir','nois','noito','tmr','noi'},'location','best',...
           'NumColumns',2,'fontsize',7);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    grid on
    hold off
    
    %%% E component
    j=2;
    isub = (i-(ifig-1)*3-1)*ncol+j;
    subplot(nrow,ncol,isub)
    box on
    % S signal
    semilogx(pdft,10*log10(pdeq1Ess(:,i)),'-','color','k','linew',1); hold on
    semilogx(pdft,10*log10(pdeq2Ess(:,i)),'-','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pdft,10*log10(pdeq4Ess(:,i)),'-','color',[0.7 0.7 0.7],'linew',1);
    % S noise before origin
    semilogx(pdft,10*log10(pdeq1Enoi(:,i)),'-.','color','k','linew',1);
    semilogx(pdft,10*log10(pdeq2Enoi(:,i)),'-.','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pdft,10*log10(pdeq4Enoi(:,i)),'-.','color',[0.7 0.7 0.7],'linew',1);
    % tremor
    semilogx(pdft,10*log10(pdtmrE(:,i)),'--','color','k','linew',1);
    % ambient noise
    semilogx(pdft,10*log10(pdnoiE(:,i)),':','color','k','linew',1);
    
    xticks([1e0, 1e1]);
    xlim([0.5, Fs/2]);
%     ylim(f.ax(isub),[-150,-20]);
%     xlabel(f.ax(isub),'Frequency (Hz)');
%     ylabel(f.ax(isub),'PSD (dB/Hz)');
    title(strcat(stnmsel(i,:),{'.E'}));
%     legend(f.ax(isub),[p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');
    grid on
    hold off


    %%% E component
    j=3;
    isub = (i-(ifig-1)*3-1)*ncol+j;
    subplot(nrow,ncol,isub)
    box on
    % S signal
    semilogx(pdft,10*log10(pdeq1Uss(:,i)),'-','color','k','linew',1); hold on
    semilogx(pdft,10*log10(pdeq2Uss(:,i)),'-','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pdft,10*log10(pdeq4Uss(:,i)),'-','color',[0.7 0.7 0.7],'linew',1);
    % S noise before origin
    semilogx(pdft,10*log10(pdeq1Unoi(:,i)),'-.','color','k','linew',1);
    semilogx(pdft,10*log10(pdeq2Unoi(:,i)),'-.','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pdft,10*log10(pdeq4Unoi(:,i)),'-.','color',[0.7 0.7 0.7],'linew',1);
    % tremor
    semilogx(pdft,10*log10(pdtmrU(:,i)),'--','color','k','linew',1);
    % ambient noise
    semilogx(pdft,10*log10(pdnoiU(:,i)),':','color','k','linew',1);
    
    xticks([1e0, 1e1]);
    xlim([0.5, Fs/2]);
%     ylim(f.ax(isub),[-150,-20]);
%     xlabel(f.ax(isub),'Frequency (Hz)');
%     ylabel(f.ax(isub),'PSD (dB/Hz)');
    title(strcat(stnmsel(i,:),{'.U'}));
%     legend(f.ax(isub),[p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');
    grid on
    hold off
    
    if mod(i,3)==0 || i== size(stnmsel,1)
        % save figure
        print(f.fig,'-dpdf',strcat(figpath,'/PSD.pd.',num2str(ifig),prefix,'.dep',num2str(depran),...
             '.pdf'));
    end
end




%% plot pmtm & pchave PSD 

ifig = 0;
for i = 1: size(stnmsel,1)     % each station
    
    if mod(i,3)==1
       
        f2.fig=figure;
        f2.fig.Renderer='Painters';
        widin = 8.3;  % maximum width allowed is 8.5 inches
        htin = 9;   % maximum height allowed is 11 inches
        set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

        ifig = ifig+1;
    end
    
    nrow = 3;     % number of stations
    ncol = 3;   % number of components
    
    %%% N component
    
    j=1;
    isub = (i-(ifig-1)*nrow-1)*ncol+j;
    subplot(nrow,ncol,isub)
    box on
    % S signal
    semilogx(pmft,10*log10(pmeq1Nss(:,i)),'-','color','k','linew',1); hold on
    semilogx(pmft,10*log10(pmeq2Nss(:,i)),'-','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pmft,10*log10(pmeq4Nss(:,i)),'-','color',[0.7 0.7 0.7],'linew',1);
    % S noise before origin
    semilogx(pmft,10*log10(pmeq1Nnoi(:,i)),'-.','color','k','linew',1);
    semilogx(pmft,10*log10(pmeq2Nnoi(:,i)),'-.','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pmft,10*log10(pmeq4Nnoi(:,i)),'-.','color',[0.7 0.7 0.7],'linew',1);
    % tremor
    semilogx(pcft,10*log10(pctmrN(:,i)),'--','color','k','linew',1);
    % ambient noise
    semilogx(pcft,10*log10(pcnoiN(:,i)),':','color','k','linew',1);
    
    xticks([1e0, 1e1]);
    xlim([0.5, Fs/2]);
%     ylim(f.ax(isub),[-150,-20]);
%     xlabel(f.ax(isub),'Frequency (Hz)');
%     ylabel(f.ax(isub),'PSD (dB/Hz)');
    title(strcat(stnmsel(i,:),{'.N pmtm & pchave'}));
    lgd = legend({'evtir','evts','evtto','noiir','nois','noito','tmr','noi'},'location','best',...
           'NumColumns',2,'fontsize',7);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    grid on
    hold off
    
    %%% E component
    j=2;
    isub = (i-(ifig-1)*3-1)*ncol+j;
    subplot(nrow,ncol,isub)
    box on
    % S signal
    semilogx(pmft,10*log10(pmeq1Ess(:,i)),'-','color','k','linew',1); hold on
    semilogx(pmft,10*log10(pmeq2Ess(:,i)),'-','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pmft,10*log10(pmeq4Ess(:,i)),'-','color',[0.7 0.7 0.7],'linew',1);
    % S noise before origin
    semilogx(pmft,10*log10(pmeq1Enoi(:,i)),'-.','color','k','linew',1);
    semilogx(pmft,10*log10(pmeq2Enoi(:,i)),'-.','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pmft,10*log10(pmeq4Enoi(:,i)),'-.','color',[0.7 0.7 0.7],'linew',1);
    % tremor
    semilogx(pcft,10*log10(pctmrE(:,i)),'--','color','k','linew',1);
    % ambient noise
    semilogx(pcft,10*log10(pcnoiE(:,i)),':','color','k','linew',1);
    
    xticks([1e0, 1e1]);
    xlim([0.5, Fs/2]);
%     ylim(f.ax(isub),[-150,-20]);
%     xlabel(f.ax(isub),'Frequency (Hz)');
%     ylabel(f.ax(isub),'PSD (dB/Hz)');
    title(strcat(stnmsel(i,:),{'.E'}));
%     legend(f.ax(isub),[p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');
    grid on
    hold off


    %%% E component
    j=3;
    isub = (i-(ifig-1)*3-1)*ncol+j;
    subplot(nrow,ncol,isub)
    box on
    % S signal
    semilogx(pmft,10*log10(pmeq1Uss(:,i)),'-','color','k','linew',1); hold on
    semilogx(pmft,10*log10(pmeq2Uss(:,i)),'-','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pmft,10*log10(pmeq4Uss(:,i)),'-','color',[0.7 0.7 0.7],'linew',1);
    % S noise before origin
    semilogx(pmft,10*log10(pmeq1Unoi(:,i)),'-.','color','k','linew',1);
    semilogx(pmft,10*log10(pmeq2Unoi(:,i)),'-.','color',[0.5 0.5 0.5],'linew',1);
    semilogx(pmft,10*log10(pmeq4Unoi(:,i)),'-.','color',[0.7 0.7 0.7],'linew',1);
    % tremor
    semilogx(pcft,10*log10(pctmrU(:,i)),'--','color','k','linew',1);
    % ambient noise
    semilogx(pcft,10*log10(pcnoiU(:,i)),':','color','k','linew',1);
    
    xticks([1e0, 1e1]);
    xlim([0.5, Fs/2]);
%     ylim(f.ax(isub),[-150,-20]);
%     xlabel(f.ax(isub),'Frequency (Hz)');
%     ylabel(f.ax(isub),'PSD (dB/Hz)');
    title(strcat(stnmsel(i,:),{'.U'}));
%     legend(f.ax(isub),[p(1) p(2) p(3)],{'PGC','SSIB','SILB'},'location','southwest');
    grid on
    hold off
    
    if mod(i,3)==0 || i== size(stnmsel,1)
        % save figure
        print(f2.fig,'-dpdf',strcat(figpath,'/PSD.pmpc.',num2str(ifig),prefix,'.dep',num2str(depran),'.pdf'));
    end
    
end


%% plot these events on map, and also the stations that I have data
f3.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f3.ax = gca;
hold(f3.ax,'on');
grid(f3.ax, 'on');
f3.ax.Box = 'on';
if regflag == 1  % means western shikoku
    plot(f3.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
         [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f3.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
axis(f3.ax, 'equal');
if regflag == 1
    xlim(f3.ax,[131 135]);
    ylim(f3.ax,[32.3 35.2]);
elseif regflag == 2
    xlim(f3.ax,[134.5 137.2]);
    ylim(f3.ax,[33 35.5]);
end
coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(f3.ax,coast(:,1),coast(:,2),'black','linew',0.5);

scatter(f3.ax,obara(:,6),obara(:,7),2,[0 1 1],'filled','o');

p1=scatter(f3.ax,evtirtest(:,8),evtirtest(:,9),40,'k','filled','o','markeredgec','k');

p2=scatter(f3.ax,evtstest(2,8),evtstest(2,9),40,[0.5 0.5 0.5],'filled','o','markeredgec','k');

p3=scatter(f3.ax,evttotest(:,8),evttotest(:,9),40,[0.7 0.7 0.7],'filled','o','markeredgec','k');


% scatter(evtsmin(:,8),evtsmin(:,9),8,'b','filled','o','markere','k');


fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationlist'), 'r');
format = '%s %f %f %f \n';
catcell = textscan(fid, format);
fclose(fid);
stnm = char(catcell{1,1});    % station name, string
stlo = catcell{1,2};    % station longitude
stla = catcell{1,3};    % station latitude
stel = catcell{1,4};    % station elevation

scatter(f3.ax,stlo,stla,60,'y','filled','^','markeredgec','k');

scatter(f3.ax,stlosel,stlasel,60,'r','filled','^','markeredgec','k');

text(f3.ax,stlo, stla+0.03, stnm);

evtirid = strcat(num2str(evtirtest(1,1)),{' '},num2str(evtirtest(1,5)),':',num2str(evtirtest(1,6)),':',...
                 num2str(evtirtest(1,7)),{', mag= '},num2str(evtirtest(1,11)));
evtsid = strcat(num2str(evtstest(2,1)),{' '},num2str(evtstest(2,5)),':',num2str(evtstest(2,6)),':',...
                num2str(evtstest(2,7)),{', mag= '},num2str(evtstest(2,11)));
evttoid = strcat(num2str(evttotest(1,1)),{' '},num2str(evttotest(1,5)),':',num2str(evttotest(1,6)),':',...
                 num2str(evttotest(1,7)),{', mag= '},num2str(evttotest(1,11)));
legend(f3.ax,[p1,p2,p3],{string(strcat({'evt within 2hr of tmr, '},evtirid)), string(evtsid), ...
       string(evttoid)},'location','southeast');


% save figure
print(f3.fig,'-dpdf',strcat(figpath,'/selsta_evt_loc.',prefix,'.dep',num2str(depran),'.pdf'));










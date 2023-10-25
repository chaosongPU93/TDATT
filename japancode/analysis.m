% function analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main script to analyze the spatio-temporal tremor and regular
% events in the Japanese region (now in kii and shikoku) 
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
% First created date:   2019/12/27
% Last modified date:   2020/02/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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

%% plot the occuring frequency(times)--magnitude relation of regular event 
evtmag = evtall(:,11);
f1.fig=figure;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f1.ax=pltevtmagfreq(evtmag,0.5);

% save figure
print(f1.fig,'-dpdf',strcat(figpath,'/evt.mag_freq.',prefix,'.dep',num2str(depran),'.pdf'));


%% plot the spatial distribution of the regular events with different magnitudes
f5.fig=figure;
f5.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8.5;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f5.ax(isub) = subplot(nrow,ncol,isub);
end
% reposition
set(f5.ax(1), 'position', [ 0.08, 0.55, 0.4, 0.4], 'unit','normalized');
set(f5.ax(2), 'position', [ 0.56, 0.55, 0.4, 0.4], 'unit','normalized');
set(f5.ax(3), 'position', [ 0.08, 0.06, 0.4, 0.4], 'unit','normalized');
set(f5.ax(4), 'position', [ 0.56, 0.06, 0.4, 0.4], 'unit','normalized');


f5 = pltevtloc(f5,regflag,evtall,obara,slab);
% save figure
print(f5.fig,'-dpdf',strcat(figpath,'/evt.spatialloc.',prefix,'.dep',num2str(depran),'.pdf'));
% print(f5.fig,'-depsc2',strcat(figpath,'/evt.spatialloc.',prefix,'.dep',num2str(depran),'.eps'));


%% obtain the tremor active dates
datetre = unique(obara(:,1));
dateact = [];   % tremor active dates
lend = length(datetre);
count = zeros(lend,1);
counttol = 0.25*24;
for i = 1: lend
    count(i) = sum(datetre(i) == obara(:,1));
    if count(i) >= counttol    % num of mins in one day
        dateact = [dateact; datetre(i)];
    end
end

%% obtain the overlapped the dates, i.e. there are regular earthquakes in those tremor active dates
%%% 1. during tremor active days 
lend = length(dateact);
evtobj = [];    % regular events during tremor active dates
for i = 1:lend
    evtobj = [evtobj; evtall(dateact(i) == evtall(:,1), :)];
end

% apply a magnitude threshold, comes from previous estimates
minmag = 1;
evtobj = evtobj(evtobj(:,11)>=minmag, :);

% also get the tremor detections during those active dates
obaraobj = [];  % tremor detections during tremor active dates
for i = 1:lend
    obaraobj = [obaraobj; obara(dateact(i) == obara(:,1),:)];
end

%%% 2. during tremor days but inactive, complement set
dateinacttre = setdiff(datetre,dateact);      % dates of no tremor
lend = length(dateinacttre);
treinacttre = [];     % tremor in inactive days
evtinacttre = [];     % regular earthquakes in inactive days
for i = 1:lend
   treinacttre = [treinacttre; obara(dateinacttre(i) == obara(:,1),:)];
   evtinacttre = [evtinacttre; evtall(dateinacttre(i) == evtall(:,1) & evtall(:,11)>=minmag, :)];
end

%%% 3. during tremor free days, corresponding is the all tremors, obara
datereg = unique(evtall(:,1));
datenotre = setdiff(datereg,datetre);
lend = length(datenotre);
evtnotre = [];     % regular earthquakes in tremor free days
for i = 1:lend
   evtnotre = [evtnotre; evtall(datenotre(i) == evtall(:,1) & evtall(:,11)>=minmag, :)];
end

%%% 4. during days that tremor is not active, corresponding is the active tremors, obaraobj
dateinact = setdiff(datereg,dateact);
lend = length(dateinact);
evtinact = [];     % regular earthquakes in tremor free days
for i = 1:lend
   evtinact = [evtinact; evtall(dateinact(i) == evtall(:,1) & evtall(:,11)>=minmag, :)];
end


%% some statistics about the results
evtobjst = sortrows(evtobj,[11,1],'descend');

dateobj = unique(evtobj(:,1));

d1 = datetime(obara(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
d2 = datetime(obara(end,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ndays = caldays(between(d1,d2,'days'))+1;

% not all the tremor active dates have regular earthquake
% if the regular ones DO NOT occur during tremor, the expectation of ratio of date would be 0
dobj2actrat = length(dateobj) / length(dateact);

% throughout the dates we examine, only a part of dates meet our needs, get the relationship between
% the dates and event ratio
% since the tremor catalog time coverage is smaller than regular ones, the actual time constaint is
% from the tremor catalog
% if lack of hf is caused by attenuation, then when there is tremor, regular is none or way less
%   ==> eobj2allrat << dobj2allrat
% if lack of hf is intrinsic, then when there is tremor, regular is as usual
%   ==> eobj2allrat ~= dobj2allrat
dobj2allrat = length(dateact) / ndays;
nevtall = sum(evtall(:,11)>=minmag); 
eobj2allrat = size(evtobj,1) / nevtall;
ntimes1 = eobj2allrat / dobj2allrat;

evtfreq = size(evtobj,1) / length(dateact);
evtfreqave = nevtall / ndays;
ntimes2 = evtfreq / evtfreqave ;


%% plot event counts against time
lend = length(dateobj);
ievtday = zeros(lend,1);
for i = 1:lend
    di = datetime(dateobj(i),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    ievtday(i) = caldays(between(d1,di,'days'))+1; 
end

nevtobj = zeros(lend,1);
for i = 1:lend
   nevtobj(i) = sum(evtobj(:,1)==dateobj(i)); 
end

f2.fig=figure;
widin = 6;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 4;
ncol = 1;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end
f2=pltevtduringtmr(f2,d1,d2,ievtday,nevtobj);
    
% save figure
print(f2.fig,'-dpdf',strcat(figpath,'/evtobj.count_time.',prefix,'.dep',num2str(depran),'.pdf'));

%% plot location of object tremor and regular event in different periods
% define and position the figure frame and axes of each plot
f3.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f3.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f3.ax(1), 'position', [ 0.08, 0.56, 0.4, 0.4]);
set(f3.ax(2), 'position', [ 0.56, 0.56, 0.4, 0.4]);
set(f3.ax(3), 'position', [ 0.08, 0.08, 0.4, 0.4]);
set(f3.ax(4), 'position', [ 0.56, 0.08, 0.4, 0.4]);
    
% subplot 1, regular events vs. tremor during tremor active days
f3.ax(1) = pltevttmrloc(f3.ax(1),regflag,slab,evtobj,obaraobj,minmag);
text(f3.ax(1),0.05,0.95,'regular events in tremor active days','unit','normalized');
text(f3.ax(1),0.05,0.9,'tremors in tremor active days','unit','normalized');

% subplot 2, regular events vs. tremor during tremor days, but inactive
f3.ax(2) = pltevttmrloc(f3.ax(2),regflag,slab,evtinacttre,treinacttre,minmag);
text(f3.ax(2),0.05,0.95,'regular events in tremor days, but inactive','unit','normalized');
text(f3.ax(2),0.05,0.9,'tremors in tremor days, but inactive','unit','normalized');

% subplot 3, regular events in tremor free days vs. all tremor
f3.ax(3) = pltevttmrloc(f3.ax(3),regflag,slab,evtnotre,obara,minmag);
text(f3.ax(3),0.05,0.95,'regular events in tremor free days','unit','normalized');
text(f3.ax(3),0.05,0.9,'all tremors','unit','normalized');

% subplot 4, regular events in tremor inactive days vs tremor in active days
f3.ax(4) = pltevttmrloc(f3.ax(4),regflag,slab,evtinact,obaraobj,minmag);
text(f3.ax(4),0.05,0.95,'regular events in tremor inactive days','unit','normalized');
text(f3.ax(4),0.05,0.9,'tremors in tremor active days','unit','normalized');

% save figure
print(f3.fig,'-dpdf',strcat(figpath,'/evttmrobj.spatialloc.',prefix,'.dep',num2str(depran),'.pdf'));


%% plot the entire tremor catalog inside that region including ETS estimate (counts perday vs. time)
% to get a sense of the major ETS periods
lend = length(datetre);
itmrday = zeros(lend,1);
for i = 1:lend
%     tmp = num2str(dateobj(i));
%     tmp2 = strcat(tmp(1:4),'-',tmp(5:6),'-',tmp(7:8)); 
    di = datetime(datetre(i),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    itmrday(i) = caldays(between(d1,di,'days'))+1; 
end

[stday, edday] = ETSestimate(regflag);
ndayets = edday-stday+1;
for i = 1: length(ndayets) 
    if ndayets(i) <7
        disp(edday(i));
    end
end

f4.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f4.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 4;
ncol = 1;
for isub = 1:nrow*ncol
    f4.ax(isub) = subplot(nrow,ncol,isub);
end
f4=pltETS(f4,d1,d2,itmrday,count,stday,edday);

% save figure
print(f4.fig,'-dpdf',strcat(figpath,'/tmr.ETStempo.',prefix,'.dep',num2str(depran),'.pdf'));






















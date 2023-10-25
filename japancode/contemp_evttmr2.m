% function contemp_evttmr2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main script to analyze the spatio-temporal evolution of tremor
% and regular events occurred in the same ETS period, to find the nearly, and/or
% co-located tremor and eqs, save those candidate events for further analysis
%
% version 2, different logic, treat events inside or outside tremor region separately
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
%       3. Kii is tested, Shikoku has NOT
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/26
% Last modified date:   2020/02/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
% close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = '/home/data2/chaosong/shikoku_kii';
tmrpath = strcat(workpath,'/etscat');
evtpath = strcat(workpath,'/jmacat02-18');
figpath = strcat(workpath,'/figs');

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

%% find the boundary defined by the tremor
% use function boundary to obtain the boundary index and area with default shrinking factor 0.5
[indbd,areain] = boundary(obara(:,6),obara(:,7),0.5);
tmrbd = [obara(indbd,6) obara(indbd,7)];

if regflag == 1  % means western shikoku
    rectx = [134.1 134.7 132.11 131.51]';
    recty = [35.1 33.7 32.59 33.99]';
elseif regflag == 2  % means kii pen
    rectx = [134.6 135.2 137.0 136.4]';
    recty = [34.1 33.2 34.4 35.3]';
end
[indrect,areaall] = boundary(rectx,recty,0);
rectbd = [rectx(indrect) recty(indrect)];

%% get regular events vertically within the depran below, horizontally within/out the tremor site
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


%% 
% imagine there is a tremor detection at a time at a location, what is the detection rate inside
% a horizontal radius within a time range relative to the tremor timing; in the same time range,
% what is the rate for the region outside that circle, but inside tremor region boundary; And what
% is the case for the region outside tremor region boundary

dttol = 2;      % max half time difference, in hr
dloctol = 0.1;  % max half location difference, in deg

% rate of events inside radius, outside radius, outside tremor boundary
% events inside radius, outside radius, outside tremor boundary
% tremors inside radius, outside radius, outside tremor boundary
[~,~,~,evtir,evtor,evtot,tmrir,tmror,tmrot] = ...
         evtratehr(evtinall,evtoutall,obara,1,1,dttol,dloctol,7);

% add a marker to reference events, in order to track them correctly
ind = 1:size(evtir,1);
evtir = [evtir ind'];
ind = 1:size(evtor,1);
evtor = [evtor ind'];
ind = 1:size(evtot,1);
evtot = [evtot ind'];
     
% set a empirical magnitude threshold from previous analysis 
minmag = 1.0;   % min magnitude
evtirmin = evtir(evtir(:,11)>=minmag, :);  
evtormin = evtor(evtor(:,11)>=minmag, :);
evtotmin = evtot(evtot(:,11)>=minmag, :);



%% find the some other events for comparison
%%% 1. that are spatially close, but in tremor free days
%%% 2. that are timely close, but outside the in radius regions
%%% 3. that are timely close, but outside the entire tremor regions

magtol = 5;   % max magnitude difference
dloctol = 0.1;  % max half location difference, in deg
deptol = 5;     % max depth difference, in km
dttol = 12;      % max half time difference, in hr

[evts,evtti,evtto] = evtouttmr(evtinall,evtoutall,evtir,obara,tmrir,magtol,dloctol,deptol,dttol);



evtsmin = evts(evts(:,11)>=minmag, :);
evttimin = evtti(evtti(:,11)>=minmag, :);
evttomin = evtto(evtto(:,11)>=minmag, :);


%% plot these events on map, and also the stations that I have data
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f.ax = gca;
hold(f.ax,'on');
grid(f.ax, 'on');
if regflag == 1  % means western shikoku
    plot(f.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
         [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
axis(f.ax, 'equal');
if regflag == 1
    xlim(f.ax,[131 135]);
    ylim(f.ax,[32.3 35.2]);
elseif regflag == 2
    xlim(f.ax,[134.5 137.2]);
    ylim(f.ax,[33 35.5]);
end
coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(f.ax,coast(:,1),coast(:,2),'black','linew',0.5);

lat = 31:0.01:36;
lon = 131:0.01:138;
[longrd, latgrd] = meshgrid(lon,lat);
depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');
contour(f.ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.8 0.8 0.8]);

scatter(f.ax,obara(:,6),obara(:,7),2,[0.3 0.3 0.3],'filled','o');

scatter(f.ax,evtall(:,8),evtall(:,9),2,[0.5 0.8 1],'filled','o');

scatter(f.ax,evtsmin(:,8),evtsmin(:,9),12,[0 0.9 1],'filled','o','markeredgec','k');

scatter(f.ax,evttimin(:,8),evttimin(:,9),15,'g','filled','o','markeredgec','k');

scatter(f.ax,evttomin(:,8),evttomin(:,9),20,[0.7 0.7 0.7],'filled','o','markeredgec','k');


scatter(f.ax,evtirmin(:,8),evtirmin(:,9),40,'r','filled','o','markeredgec','k');


% scatter(evtsmin(:,8),evtsmin(:,9),8,'b','filled','o','markere','k');


fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationlist'), 'r');
format = '%s %f %f %f \n';
catcell = textscan(fid, format);
fclose(fid);
stnm = char(catcell{1,1});    % station name, string
stlo = catcell{1,2};    % station longitude
stla = catcell{1,3};    % station latitude
stel = catcell{1,4};    % station elevation

scatter(f.ax,stlo,stla,60,'y','filled','^','markeredgec','k');

text(f.ax,stlo, stla+0.03, stnm);


%% for test purpose, choose some stations
ind = [6;10;13;21;23;24;30;31];
stnmsel = stnm(ind,:);
stlosel = stlo(ind);
stlasel = stla(ind);
stelsel = stel(ind);

scatter(f.ax,stlosel,stlasel,60,'k','filled','^','markeredgec','k');


fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationselnm'),'w+');
for i = 1: size(stnmsel,1)
    fprintf(fid,'%s \n',stnmsel(i,:));
end
fclose(fid);

fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationlistsel'),'w+');
for i = 1: size(stnmsel,1)
    fprintf(fid,'%s %f %f %f \n',stnmsel(i,:),stlosel(i),stlasel(i),stelsel(i));
end
fclose(fid);

% all the dates
dtevtirmin = unique(evtirmin(:,1));
dtevttomin = unique(evttomin(:,1));
dtevtsmin = unique(evtsmin(:,1));
dtevttimin = unique(evttimin(:,1));
dateall = [dtevtirmin;dtevttomin;dtevtsmin;dtevttimin ];

fid = fopen(strcat('/home/data2/chaosong/japan_chao/selectdate'),'w+');
fprintf(fid,'%d \n',dateall');
fclose(fid);


%% save the events data for future use 
save(strcat(evtpath,'/objevts_forspectra.',prefix,'dep',num2str(depran),'.mat'),'evtir','evts',...
     'evtti','evtto');




















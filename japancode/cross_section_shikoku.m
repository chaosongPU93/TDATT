% function cross_section_shikoku.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to plot the cross-section of the earthquake distribution for shikoku region
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
% First created date:   2020/03/24
% Last modified date:   2020/03/24
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
sacpath = '/home/data2/chaosong/japan_chao';

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

% event type, 1 for natural EQ, 5 for LFE
evttype = zeros(size(evtall,1),1);
for i = 1: size(evtall,1)
    evttype(i) = str2double(event(i).evttp);
end
lfe = evtall(evttype==5, :);    % lfe events
perclfe = size(lfe,1)/size(evtall,1);

% magnitude type, V-velocity mag; v-as per V, but for 2 or 3 stas; 
mag1type = char(size(evtall,1),1);
for i = 1: size(evtall,1)
    mag1type(i,:) = event(i).mag1tp;
end
evtmagur = evtall(mag1type=='v', :);   % events whose mag is slightly unreliable
percmagur = size(evtmagur,1)/size(evtall,1);

% hypo location precison flag, 1 for depth-free method, 7 for poor solution (before 2016)
hypprec = char(size(evtall,1),1);
for i = 1: size(evtall,1)
    hypprec(i) = event(i).hyplocprec;
end
evthyppoor = evtall(hypprec=='7', :);
perchyppoor = size(evthyppoor,1)/size(evtall,1);

% num of stas used for hypocenter determination
numsta = zeros(size(evtall,1),1);
for i = 1: size(evtall,1)
    numsta(i) = event(i).nsta;
end
minnumsta = min(numsta);    % see how small can the number of stations be

% Hypocenter determination flag, K for high precision; S/s for low/inferior precision
hypdet = char(size(evtall,1),1);
for i = 1: size(evtall,1)
    hypdet(i) = event(i).hypdetflag;
end 
evthypur = evtall(hypdet=='S' | hypdet=='s', :);   % events whose mag is slightly unreliable
perchypur = size(evthypur,1)/size(evtall,1);


%% check the event catalog more carefully, retain the natural EQ only
evtall = evtall(evttype==1, :);     % 1 is natural EQ



%% plot these events on map, and also the stations that I have data
f.fig=figure;
f.fig.Renderer='Painters';
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

scatter(f.ax,evtall(:,8),evtall(:,9),2,[0.5 0.8 1],'filled','o');

scatter(f.ax,obara(:,6),obara(:,7),2,[0.3 0.3 0.3],'filled','o');

fid = fopen(strcat(workpath,'/NIED_SeismicStation_20200324.csv'), 'r');
format = '%s %s %s %s %s %f %f %f %f %f %s %s \n';
catcell = textscan(fid, format, 'delimiter',',','headerlines',1);
fclose(fid);
stnm = char(catcell{1,3});    % station name, string
stlo = catcell{1,7};    % station longitude
stla = catcell{1,6};    % station latitude
stel = catcell{1,8};    % station elevation

%%% bound type 0, the boundary of tremor
if regflag == 1  % means western shikoku
    bnd0 = [134.1 35.1;
            134.7 33.7;
            132.11 32.59;
            131.51 33.99;
            134.1 35.1];
elseif regflag == 2  % means kii pen
    bnd0 = [134.6 34.1;
            135.2 33.2;
            137 34.4;
            136.4 35.3;
            134.6 34.1];
end

% only show the station inside bound 0
[is,ion] = inpolygon(stlo,stla,bnd0(:,1),bnd0(:,2));
isinbnd0 = is | ion;
stnm = stnm(isinbnd0 == 1, :);
stlo = stlo(isinbnd0 == 1);
stla = stla(isinbnd0 == 1);
stel = stel(isinbnd0 == 1);

scatter(f.ax,stlo, stla,60,'y','filled','^','markeredgec','k');

text(f.ax,stlo, stla+0.03, stnm, 'fontsize',6);


%% set the further boundary to regular events 
% from the distribution of the regular events, it seems proper to separate the whole region into 2
% parts, both of them are outlined by hand

%%% bound 1, the NE part
bnd1 = [134.3197   34.1143;
        134.5787   33.9625;
        134.6377   33.7597;
        134.3440   33.5630;
        134.0303   33.4262;
        133.8752   33.5045;
        133.5283   33.9357;
        133.6370   34.0510;
        133.8827   34.1191;
        134.3197   34.1143];

%%% bound type 2, irregular shape hand picked by choosing the co-located region as tremors
bnd2 = [132.1470   33.4843;
        132.5841   33.6345;
        132.8415   33.4803;
        132.8503   33.2645;
        132.75     32.87;
        132.11     32.59;
        131.9102   33.0760;
        131.8575   33.3213;
        132.1470   33.4843];

% events & tremors inside bound 1
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),bnd1(:,1),bnd1(:,2));
isinbnd1 = is | ion;
evtbnd1 = evtall(isinbnd1 == 1, :);
[is,ion] = inpolygon(obara(:,6),obara(:,7),bnd1(:,1),bnd1(:,2));
isinbnd1 = is | ion;
tmrbnd1 = obara(isinbnd1 == 1, :);

% events & tremors inside bound 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),bnd2(:,1),bnd2(:,2));
isinbnd2 = is | ion;
evtbnd2 = evtall(isinbnd2 == 1, :);
[is,ion] = inpolygon(obara(:,6),obara(:,7),bnd2(:,1),bnd2(:,2));
isinbnd2 = is | ion;
tmrbnd2 = obara(isinbnd2 == 1, :);

plot(f.ax,bnd1(:,1),bnd1(:,2),'color','b','linew',2);
plot(f.ax,bnd2(:,1),bnd2(:,2),'color','b','linew',2);


%% cross-scetion profile to cluster 1 & 2
% csp1x = linspace(134.1,133.7,1000);
% csp1y = linspace(33.6,34.9,1000);
% csp1 = [csp1x' csp1y'];
% 
% csp2x = linspace(132.6,131.7,1000);
% csp2y = linspace(33.0,33.8,1000);
% csp2 = [csp2x' csp2y'];

% cross-section region 1 (rect)
csp1 = [134 33.35;
        133.6 34.25;
        133.5178 34.2134;
        133.9178 33.3134;
        134 33.35];

% the middle line of the cross-section
csp1midx = linspace(0.5*(csp1(4,1)+csp1(5,1)),0.5*(csp1(2,1)+csp1(3,1)),1000);
csp1midy = linspace(0.5*(csp1(4,2)+csp1(5,2)),0.5*(csp1(2,2)+csp1(3,2)),1000);
csp1mid = [csp1midx' csp1midy'];

% cross-section region 2 (rect)
csp2 = [132.75 32.87;
        131.7 33.8;
        131.6402 33.7327;
        132.6902 32.8027;
        132.75 32.87];
    
% the middle line of the cross-section 2
csp2midx = linspace(0.5*(csp2(4,1)+csp2(5,1)),0.5*(csp2(2,1)+csp2(3,1)),1000);
csp2midy = linspace(0.5*(csp2(4,2)+csp2(5,2)),0.5*(csp2(2,2)+csp2(3,2)),1000);
csp2mid = [csp2midx' csp2midy'];

% cross-section region 2-2 (rect)
csp3 = [133.1 33.25;
        132.05 34.16;
        131.9171 34.0105;
        132.9671 33.1005;
        133.1 33.25];

% the middle line of the cross-section 2-2
csp3midx = linspace(0.5*(csp3(4,1)+csp3(5,1)),0.5*(csp3(2,1)+csp3(3,1)),1000);
csp3midy = linspace(0.5*(csp3(4,2)+csp3(5,2)),0.5*(csp3(2,2)+csp3(3,2)),1000);
csp3mid = [csp3midx' csp3midy'];

% cross-section region 1-2 (rect)
csp4 = [133.95 33.5;
        133.95 34.3;
        133.86 34.3;
        133.86 33.5;
        133.95 33.5];    

% the middle line of the cross-section 1-2
csp4midx = linspace(0.5*(csp4(4,1)+csp4(5,1)),0.5*(csp4(2,1)+csp4(3,1)),1000);
csp4midy = linspace(0.5*(csp4(4,2)+csp4(5,2)),0.5*(csp4(2,2)+csp4(3,2)),1000);
csp4mid = [csp4midx' csp4midy'];

% assume tremors are on the slab interface, interpolate to get their depth
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
obara(:,8) = F(obara(:,6),obara(:,7));

patch(f.ax,csp1(:,1),csp1(:,2),'r','Facealpha',0.4);
plot(f.ax,csp1mid(:,1),csp1mid(:,2),'r','linew',1.5);
patch(f.ax,csp2(:,1),csp2(:,2),'r','Facealpha',0.4);
plot(f.ax,csp2mid(:,1),csp2mid(:,2),'r','linew',1.5);
% patch(f.ax,csp3(:,1),csp3(:,2),'r','Facealpha',0.4);
% plot(f.ax,csp3mid(:,1),csp3mid(:,2),'r','linew',1.5);
patch(f.ax,csp4(:,1),csp4(:,2),'r','Facealpha',0.4);
plot(f.ax,csp4mid(:,1),csp4mid(:,2),'r','linew',1.5);
text(f.ax,133.5,34.3,'A','FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,133.95,33.25,"A'",'FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,131.55,33.85,'B','FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,132.75,32.75,"B'",'FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,131.8,34.1,'C','FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,133.05,33.10,"C'",'FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,133.95,34.4,'D','FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,134.0,33.4,"D'",'FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
title(f.ax,prefix);
f.ax.Box='on';

% save figure
print(f.fig,'-dpdf',strcat(figpath,'/evt_tmr_sta_bounds_csection.',prefix,'.dep',...
      num2str(depran),'.pdf'));

  
%% project all events & tremor onto it  
%%%% for bound 1 events and their cross-section profile
% get events & tremors inside the profile, regardless of bound1 and 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),csp1(:,1),csp1(:,2));
isincsp = is | ion;
evtcsp1 = evtall(isincsp == 1, :);

[is,ion] = inpolygon(obara(:,6),obara(:,7),csp1(:,1),csp1(:,2));
isincsp = is | ion;
tmrcsp1 = obara(isincsp == 1, :);

% define the origin as the lower right point
lon0 = csp1mid(1,1);
lat0 = csp1mid(1,2);

% convert the cross-section profile to relative distance 
[csp1midr(:,1),csp1midr(:,2)] = absloc2relaloc(csp1mid(:,1),csp1mid(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(csp1midr(end,2),csp1midr(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
csp1midr(:,3) = F(csp1mid(:,1),csp1mid(:,2));
% rotate the coordinates to profile axis 
[csp1midr(:,4),csp1midr(:,5)] = coordinate_rot(csp1midr(:,1),csp1midr(:,2),rotang,0,0);

% convert events' locations to relative distance
[dx,dy] = absloc2relaloc(evtcsp1(:,8),evtcsp1(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtcsp1(:,13),evtcsp1(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(tmrcsp1(:,6),tmrcsp1(:,7),lon0,lat0);
[tmrcsp1(:,9),tmrcsp1(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stx1r,sty1r] = coordinate_rot(dx,dy,rotang,0,0);
%%%%

%%%% for bound 2 events and their cross-section profile
% get events & tremors inside the profile, regardless of bound1 and 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),csp2(:,1),csp2(:,2));
isincsp = is | ion;
evtcsp2 = evtall(isincsp == 1, :);

[is,ion] = inpolygon(obara(:,6),obara(:,7),csp2(:,1),csp2(:,2));
isincsp = is | ion;
tmrcsp2 = obara(isincsp == 1, :);

% define the origin as the lower right point
lon0 = csp2mid(1,1);
lat0 = csp2mid(1,2);

% convert the cross-section profile to relative distance 
[csp2midr(:,1),csp2midr(:,2)] = absloc2relaloc(csp2mid(:,1),csp2mid(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(csp2midr(end,2),csp2midr(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
csp2midr(:,3) = F(csp2mid(:,1),csp2mid(:,2));
% rotate the coordinates to profile axis 
[csp2midr(:,4),csp2midr(:,5)] = coordinate_rot(csp2midr(:,1),csp2midr(:,2),rotang,0,0);

% convert events' locations to relative distance
[dx,dy] = absloc2relaloc(evtcsp2(:,8),evtcsp2(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtcsp2(:,13),evtcsp2(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(tmrcsp2(:,6),tmrcsp2(:,7),lon0,lat0);
[tmrcsp2(:,9),tmrcsp2(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stx2r,sty2r] = coordinate_rot(dx,dy,rotang,0,0);
%%%%


%%%% for events in csp 3 
% get events & tremors inside the profile, regardless of bound1 and 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),csp3(:,1),csp3(:,2));
isincsp = is | ion;
evtcsp3 = evtall(isincsp == 1, :);

[is,ion] = inpolygon(obara(:,6),obara(:,7),csp3(:,1),csp3(:,2));
isincsp = is | ion;
tmrcsp3 = obara(isincsp == 1, :);

% define the origin as the lower right point
lon0 = csp3mid(1,1);
lat0 = csp3mid(1,2);

% convert the cross-section profile to relative distance 
[csp3midr(:,1),csp3midr(:,2)] = absloc2relaloc(csp3mid(:,1),csp3mid(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(csp3midr(end,2),csp3midr(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
csp3midr(:,3) = F(csp3mid(:,1),csp3mid(:,2));
% rotate the coordinates to profile axis 
[csp3midr(:,4),csp3midr(:,5)] = coordinate_rot(csp3midr(:,1),csp3midr(:,2),rotang,0,0);

% convert events' locations to relative distance
[dx,dy] = absloc2relaloc(evtcsp3(:,8),evtcsp3(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtcsp3(:,13),evtcsp3(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(tmrcsp3(:,6),tmrcsp3(:,7),lon0,lat0);
[tmrcsp3(:,9),tmrcsp3(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stx3r,sty3r] = coordinate_rot(dx,dy,rotang,0,0);
%%%%


%%%% for events in csp 4 
% get events & tremors inside the profile, regardless of bound1 and 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),csp4(:,1),csp4(:,2));
isincsp = is | ion;
evtcsp4 = evtall(isincsp == 1, :);

[is,ion] = inpolygon(obara(:,6),obara(:,7),csp4(:,1),csp4(:,2));
isincsp = is | ion;
tmrcsp4 = obara(isincsp == 1, :);

% define the origin as the lower right point
lon0 = csp4mid(1,1);
lat0 = csp4mid(1,2);

% convert the cross-section profile to relative distance 
[csp4midr(:,1),csp4midr(:,2)] = absloc2relaloc(csp4mid(:,1),csp4mid(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(csp4midr(end,2),csp4midr(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
csp4midr(:,3) = F(csp4mid(:,1),csp4mid(:,2));
% rotate the coordinates to profile axis 
[csp4midr(:,4),csp4midr(:,5)] = coordinate_rot(csp4midr(:,1),csp4midr(:,2),rotang,0,0);

% convert events' locations to relative distance
[dx,dy] = absloc2relaloc(evtcsp4(:,8),evtcsp4(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtcsp4(:,13),evtcsp4(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(tmrcsp4(:,6),tmrcsp4(:,7),lon0,lat0);
[tmrcsp4(:,9),tmrcsp4(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stx4r,sty4r] = coordinate_rot(dx,dy,rotang,0,0);
%%%%


%% plot the cross-section profile
% x axis is along the profile, y axis is depth

f2.fig=figure;
widin = 10;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

% plot events in bound1
f2.ax(1) = subplot(221);
ax = f2.ax(1);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp1midr(:,5),csp1midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp1(:,10),tmrcsp1(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp1(evtcsp1(:,11)<1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,[0.5 0.8 1],'filled','o');
scatter(ax,sty1r,zeros(length(sty1r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty1r,5*ones(length(sty1r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,50]);
xlim(ax,[0,100]);
text(ax,60,15,"evts inside A-A', mag<1");
hold(ax,'off');


f2.ax(2) = subplot(222);
ax = f2.ax(2);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp2midr(:,5),csp2midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp2(:,10),tmrcsp2(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp2(evtcsp2(:,11)<1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,[0.5 0.8 1],'filled','o');
scatter(ax,sty2r,zeros(length(sty2r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty2r,5*ones(length(sty2r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,90]);
xlim(ax,[0,140]);
text(ax,80,15,"evts inside B-B'', mag<1");
hold(ax,'off');


f2.ax(3) = subplot(223);
ax = f2.ax(3);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp1midr(:,5),csp1midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp1(:,10),tmrcsp1(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp1(evtcsp1(:,11)>=1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,'b','filled','o');
scatter(ax,sty1r,zeros(length(sty1r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty1r,5*ones(length(sty1r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,50]);
xlim(ax,[0,100]);
text(ax,60,15,"evts inside A-A', mag>=1");
hold(ax,'off');


f2.ax(4) = subplot(224);
ax = f2.ax(4);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp2midr(:,5),csp2midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp2(:,10),tmrcsp2(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp2(evtcsp2(:,11)>=1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,'b','filled','o');
scatter(ax,sty2r,zeros(length(sty2r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty2r,5*ones(length(sty2r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,90]);
xlim(ax,[0,140]);
text(ax,80,15,"evts inside B-B'', mag>=1");
hold(ax,'off');

supertit(f2.ax,prefix);

% save figure
print(f2.fig,'-depsc2',strcat(figpath,'/evt_tmr_csection_view12.',prefix,'.dep',...
      num2str(depran),'.eps'));


%%%%%% for csp 3 and 4
f2.fig=figure;
widin = 10;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

% plot events in bound1
f2.ax(1) = subplot(221);
ax = f2.ax(1);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp3midr(:,5),csp3midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp3(:,10),tmrcsp3(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp3(evtcsp3(:,11)<1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,[0.5 0.8 1],'filled','o');
scatter(ax,sty3r,zeros(length(sty3r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty3r,5*ones(length(sty3r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,50]);
xlim(ax,[0,100]);
text(ax,60,15,"evts inside C-C', mag<1");
hold(ax,'off');


f2.ax(2) = subplot(222);
ax = f2.ax(2);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp4midr(:,5),csp4midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp4(:,10),tmrcsp4(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp4(evtcsp4(:,11)<1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,[0.5 0.8 1],'filled','o');
scatter(ax,sty4r,zeros(length(sty4r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty4r,5*ones(length(sty4r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,90]);
xlim(ax,[0,140]);
text(ax,80,15,"evts inside D-D'', mag<1");
hold(ax,'off');


f2.ax(3) = subplot(223);
ax = f2.ax(3);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp3midr(:,5),csp3midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp3(:,10),tmrcsp3(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp3(evtcsp3(:,11)>=1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,'b','filled','o');
scatter(ax,sty3r,zeros(length(sty3r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty3r,5*ones(length(sty3r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,50]);
xlim(ax,[0,100]);
text(ax,60,15,"evts inside C-C', mag>=1");
hold(ax,'off');


f2.ax(4) = subplot(224);
ax = f2.ax(4);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
% axis(ax, 'equal');
plot(ax,csp4midr(:,5),csp4midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp4(:,10),tmrcsp4(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp4(evtcsp4(:,11)>=1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,'b','filled','o');
scatter(ax,sty4r,zeros(length(sty4r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty4r,5*ones(length(sty4r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,90]);
xlim(ax,[0,140]);
text(ax,80,15,"evts inside D-D'', mag>=1");
hold(ax,'off');

supertit(f2.ax,prefix);

% save figure
print(f2.fig,'-depsc2',strcat(figpath,'/evt_tmr_csection_view34.',prefix,'.dep',...
      num2str(depran),'.eps'));
  
  
%% create a 3D plot
f3.fig=figure;
widin = 12;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

% plot events in bound1
f3.ax(1) = subplot(221);
ax = f3.ax(1);
hold on
box on
grid on
color(:,:,1) = 0.1*ones(size(latgrd));
color(:,:,2) = 0.6*ones(size(latgrd));
color(:,:,3) = 0.6*ones(size(latgrd));

% surf(longrd,latgrd,depgrd,color,'facealpha',0.5); hold on
scatter3(ax,obara(:,6),obara(:,7),obara(:,8),2,[0.3 0.3 0.3],'filled','o');
tmp = evtbnd1(evtbnd1(:,11)<1, :);
scatter3(ax,tmp(:,8),tmp(:,9),tmp(:,10),2,[0.5 0.8 1],'filled','o');
xlabel(ax,'Lon (deg)');
ylabel(ax,'Lat (deg)');
zlabel(ax,'Depth (km)');
set(ax, 'ZDir','reverse');
view(ax,285,15)
hold off

f3.ax(2) = subplot(222);
ax = f3.ax(2);
hold on
box on
grid on
scatter3(ax,obara(:,6),obara(:,7),obara(:,8),2,[0.3 0.3 0.3],'filled','o');
tmp = evtbnd2(evtbnd2(:,11)<1, :);
scatter3(ax,tmp(:,8),tmp(:,9),tmp(:,10),2,[0.5 0.8 1],'filled','o');
xlabel(ax,'Lon (deg)');
ylabel(ax,'Lat (deg)');
zlabel(ax,'Depth (km)');
set(ax, 'ZDir','reverse');
view(ax,285,15)
hold off

f3.ax(3) = subplot(223);
ax = f3.ax(3);
hold on
box on
grid on
scatter3(ax,obara(:,6),obara(:,7),obara(:,8),2,[0.3 0.3 0.3],'filled','o');
tmp = evtbnd1(evtbnd1(:,11)>=1, :);
scatter3(ax,tmp(:,8),tmp(:,9),tmp(:,10),2,[0.5 0.8 1],'filled','o');
xlabel(ax,'Lon (deg)');
ylabel(ax,'Lat (deg)');
zlabel(ax,'Depth (km)');
set(ax, 'ZDir','reverse');
view(ax,285,15)
hold off

f3.ax(4) = subplot(224);
ax = f3.ax(4);
hold on
box on
grid on
scatter3(ax,obara(:,6),obara(:,7),obara(:,8),2,[0.3 0.3 0.3],'filled','o');
tmp = evtbnd2(evtbnd2(:,11)>=1, :);
scatter3(ax,tmp(:,8),tmp(:,9),tmp(:,10),2,[0.5 0.8 1],'filled','o');
xlabel(ax,'Lon (deg)');
ylabel(ax,'Lat (deg)');
zlabel(ax,'Depth (km)');
set(ax, 'ZDir','reverse');
view(ax,285,15)
hold off

supertit(f3.ax,prefix);

% save figure
savefig(f3.fig,strcat(figpath,'/evt_tmr_3dview.',prefix,'.dep',...
      num2str(depran),'.fig'));










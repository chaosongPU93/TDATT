% function evt_piercing_point_shikoku.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to obtain the piercing point of the ray path from the 
% event to station on the slab interface. First, creat a kd-tree model 
% (KDtreeSearcher) based on the slab grid data (after interpolation). Second,
% use k-nearest neighbours (knn) searcher (knnsearch or a radius search 
% rangesearch) to find k-nearest points on the slab that are close to the 
% interpolated ray path. 
% 
%
% kdtree:
% 1.Kd-trees divide your data into nodes with at most BucketSize (default
%   is 50) points per node, based on coordinates (as opposed to categories).
%
%
%
% 
% knn search:
% 1.Given a set X of n points and a distance function, k-nearest neighbor (kNN)
%   search lets you find the k closest points in X to a query point or set of 
%   points Y.
% 2.In contrast, for a positive real value r, rangesearch finds all points in X
%   that are within a distance r of each point in Y. This fixed-radius search is
%   closely related to kNN search, as it supports the same distance metrics and 
%   search classes, and uses the same search algorithms.
% 3.knnsearch does the following: Determines the node to which the query point 
%   belongs. Finds the closest k points within that node and its distance to the
%   query point. Chooses all other nodes having any area that is within the same
%   distance, in any direction, from the query point to the kth closest point. 
%   Searches nodes within that range for any points closer to the query point.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/13
% Last modified date:   2020/04/13
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
datapath = strcat(workpath,'/matsave');

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

% %%% check the event catalog more carefully, retain the natural EQ only
% evtall = evtall(evttype==1, :);     % 1 is natural EQ


% load all Hinet stations, though i don't have any data
fid = fopen(strcat(workpath,'/NIED_SeismicStation_20200324.csv'), 'r');
format = '%s %s %s %s %s %f %f %f %f %f %s %s \n';
catcell = textscan(fid, format, 'delimiter',',','headerlines',1);
fclose(fid);
stnm = char(catcell{1,3});    % station name, string
stlo = catcell{1,7};    % station longitude
stla = catcell{1,6};    % station latitude
stel = catcell{1,8};    % station elevation


%% set the further boundary to regular events 
% from the distribution of the regular events, it seems proper to separate the whole region into 2
% parts, both of them are outlined by hand

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


%% cross-scetion profile to cluster 1, project all events & tremor onto it
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

% assume tremors are on the slab interface, interpolate to get their depth
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
obara(:,8) = F(obara(:,6),obara(:,7));

%%%% for bound 1 events and their cross-section profile
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
[dx,dy] = absloc2relaloc(evtbnd1(:,8),evtbnd1(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtbnd1(:,13),evtbnd1(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(obara(:,6),obara(:,7),lon0,lat0);
[obara(:,9),obara(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stx1r,sty1r] = coordinate_rot(dx,dy,rotang,0,0);
%%%%

%%%% for bound 2 events and their cross-section profile
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
[dx,dy] = absloc2relaloc(evtbnd2(:,8),evtbnd2(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtbnd2(:,13),evtbnd2(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(obara(:,6),obara(:,7),lon0,lat0);
[obara(:,11),obara(:,12)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stx2r,sty2r] = coordinate_rot(dx,dy,rotang,0,0);
%%%%


%% keep only staions that can sample the portion above the tremor, rather than updip slab
%   Find the upper depth limit of tremors, any stations whose along profile distance
%	smaller than that at this depth limit can be discarded, because those stations only sample
%	the updip portion of the slab, which is definitely not the zone we are seeking.
[mindep, ind] = min(obara(:,8));

% selected stations for bound 1 events and tremors
indsel1 = find(sty1r>=obara(ind,10));
stasel1 = [];
for i = 1: length(indsel1)
    stasel1(i).nm = stnm(indsel1(i),:);
    stasel1(i).lo = stlo(indsel1(i));    % station longitude
    stasel1(i).la = stla(indsel1(i));    % station latitude
    stasel1(i).el = stel(indsel1(i));    % station elevation
end


% selected stations for bound 2 events and tremors
indsel2 = find(sty2r>=obara(ind,12));
stasel2 = [];
for i = 1: length(indsel2)
    stasel2(i).nm = stnm(indsel2(i),:);
    stasel2(i).lo = stlo(indsel2(i));    % station longitude
    stasel2(i).la = stla(indsel2(i));    % station latitude
    stasel2(i).el = stel(indsel2(i));    % station elevation
end


%% for boundary 1 events, calculate the piercing points of events on the slab interface 
recalflag = 1;
dslab = 0.002;
evtsel = evtbnd1;
stasel = stasel1;
kdense = 10;
% velmod = 'jma2001';
% velmod = 'ukawa83';
velmod = 'ak135';
rayinfo1 = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
ttsevt1 =rayinfo1.tts;
raypevt1 = rayinfo1.rayp;
distevt1 = rayinfo1.dist;
pierptevt1 = rayinfo1.pierpt ;
poorresevt1 = rayinfo1.poorres;
edistevt1 = rayinfo1.edist;
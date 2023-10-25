% function same_raypath_kii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to address how to find two events that share the similar ray path.
%
%   1. situation is the same for one tremor and one event, because our goal is to find
%       the tremor and event pairs that both travel through the same region above the slab
%       (assume to have high pore pressure, high vp-vs ratio, partial melting, high
%       attenuation, etc).
%   2. It might be useful to find those pairs that happen in the same times and different
%       times, (it is possible to subtract travel times to agjust the origin time difference
%       to make sure that 'they pass through one region in the same time, rather than simply
%       having the same origin time')
%   3. The most promising tool for travel path is taup, the key is to quantify the 'same path'
%       mathmetically, although plotting is also necessary. And among all tools along with taup,
%       taup_time provides the most quantitative deteails of the ray paths. As the location of
%       of source (evt and tmr catalog) and receiver, preferred phase (S or s in our case) 
%       are given, the output is the distance [deg], source depth [km], travel time [s], ray
%       parameter p= r*sin(theta)/v, constant for the same ray [s/deg], takeoff angle i, 0~180
%       [deg], angle between take off and downward direction, incident angle, angle of arrival
%       between ray and upward direction [deg].
%
% How to quantify the 'same path' (constraints)
%   1. Since the ray p is constant along the ray, all possible sources to the same staion that
%       have the same ray p are on the same path. Their hypocenters are NOT necessarily the same, 
%       thus the distance, travel time.
%   2. Since the 1D velocity model is laterally homogeneous, it means that a source even does not
%       share a same path could have the same ray p, so all possible sources must have the same
%       azimuth.
% About the timing:
%   3. The first step is to find the same path; then we need to care about the timing. As mentioned
%       above, calibrate thee travel time difference to origin time, to make sure they pass through
%       the same region at the same time
% About the stations:
%   4. Not all stations are useful. From the distribution of tremors and regular events, it seems
%       that only stations on the downdip direction of the slab could allow the paths of tremor and
%       regular events to pass through the same region above the tremor. Otherwise, the path would
%       cross the oceanic lithosphere, which is not what we want. This could be filtered through
%       looking at the along slab-perpendicular profile distribution of tremors, events and
%       stations. Find the upper depth limit of tremors, any stations whose along profile distance
%       smaller than that at this depth limit can be discarded, because those stations only sample
%       the updip portion of the slab, which is definitely not the zone we are seeking. 
%       
%
%
% Possible difficulty
%   1. The slab model used is what we can do best, as far as I know.
%   2. The velocity model is essential to the taup calculation, so a approprite local velocity model
%       seems to be very important to us. From test, the global model 'ak135' provided by taup gives
%       the best travel time prediction. So now i using that as a start, but if the method works, i
%       need to build a local model in taup.
%
%
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/03/18
% Last modified date:   2020/03/18
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
stnmsel1 = stnm(indsel1,:);    % station name, string
stlosel1 = stlo(indsel1);    % station longitude
stlasel1 = stla(indsel1);    % station latitude
stelsel1 = stel(indsel1);    % station elevation

% selected stations for bound 2 events and tremors
indsel2 = find(sty2r>=obara(ind,12));
stnmsel2 = stnm(indsel2,:);    % station name, string
stlosel2 = stlo(indsel2);    % station longitude
stlasel2 = stla(indsel2);    % station latitude
stelsel2 = stel(indsel2);    % station elevation


%% massive calculation
% time for travel time and ray p for all stations and tremors in bound 0 take 1610 s, ~0.5 h
% time for travel time and ray p for all stations and events in bound 1 take 3574 s, ~1 h

% several available velocity models
jma2001 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.taup';
ukawa83 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ukawa83.taup';

% choose one
velmod = jma2001;
% velmod = ukawa83;
% velmod = 'ak135';

% get the char velflag for saving files
words = getwords(velmod,'/');
tmp = getwords(words{end},'.');
velflag = tmp{1};


%%%%% for bound 1 events and bound 0 tremors
tmrsel = obara;
nsta1 = size(stnmsel1,1);
ntmr1 = size(tmrsel,1);
nevt1 = size(evtbnd1,1);

% recalculation flag
recalflag = 0;

if recalflag
    % calculate the ray parameter between all stations and tremors
    ttstmr1 = zeros(ntmr1, nsta1);
    rayptmr1 = zeros(ntmr1, nsta1);
    aztmr1 = zeros(ntmr1, nsta1);
    tic
    for i = 1: nsta1
        for j = 1: ntmr1
            % travel time and ray p
            tt=tauptime('mod',velmod,'dep',tmrsel(j,8),'ph','s,S','evt',[tmrsel(j,7) tmrsel(j,6)],...
                        'sta',[stlasel1(i) stlosel1(i)]);
            ttstmr1(j,i) = tt(1).time;   % choose the first S arrival
            rayptmr1(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            
        end
        
        % azimuth from station to all tremors
        [~,az] = distance(stlasel1(i),stlosel1(i),tmrsel(:,7),tmrsel(:,6));
        aztmr1(:,i) = az;
    end
    toc
    save(strcat(datapath,'/ttsraypazi_tmr_bnd1_',prefix,'_',velflag,'.mat'),'ttstmr1','rayptmr1',...
         'aztmr1');
    
    % calculate the ray parameter between all stations and events
    ttsevt1 = zeros(nevt1, nsta1);
    raypevt1 = zeros(nevt1, nsta1);
    azevt1 = zeros(nevt1, nsta1);
    tic
    for i = 1: nsta1
        for j = 1: nevt1
            % travel time and ray p
            tt=tauptime('mod',velmod,'dep',evtbnd1(j,10),'ph','s,S','evt',...
                        [evtbnd1(j,9) evtbnd1(j,8)],'sta',[stlasel1(i) stlosel1(i)]);
            ttsevt1(j,i) = tt(1).time;   % choose the first S arrival
            raypevt1(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            
        end
        
        % azimuth from station to all tremors
        [~,az] = distance(stlasel1(i),stlosel1(i),evtbnd1(:,9),evtbnd1(:,8));
        azevt1(:,i) = az;
        
    end
    toc
    save(strcat(datapath,'/ttsraypazi_evt_bnd1_',prefix,'_',velflag,'.mat'),'ttsevt1','raypevt1',...
         'azevt1');
    
else
    data = load(strcat(datapath,'/ttsraypazi_tmr_bnd1_',prefix,'_',velflag,'.mat'));
    ttstmr1 = data.ttstmr1;
    rayptmr1 = data.rayptmr1;
    aztmr1 = data.aztmr1;
    
    data = load(strcat(datapath,'/ttsraypazi_evt_bnd1_',prefix,'_',velflag,'.mat'));
    ttsevt1 = data.ttsevt1;
    raypevt1 = data.raypevt1;
    azevt1 = data.azevt1;

end
%%%%%


%%%%% for bound 2 events and bound 0 tremors
tmrsel = obara;
nsta2 = size(stnmsel2,1);
ntmr2 = size(tmrsel,1);
nevt2 = size(evtbnd2,1);

recalflag = 0;

if recalflag
    % calculate the ray parameter between all stations and tremors
    ttstmr2 = zeros(ntmr2, nsta2);
    rayptmr2 = zeros(ntmr2, nsta2);
    aztmr2 = zeros(ntmr2, nsta2);
    tic
    for i = 1: nsta2
        for j = 1: ntmr2
            % travel time and ray p
            tt=tauptime('mod',velmod,'dep',tmrsel(j,8),'ph','s,S','evt',[tmrsel(j,7) tmrsel(j,6)],...
                        'sta',[stlasel2(i) stlosel2(i)]);
            ttstmr2(j,i) = tt(1).time;   % choose the first S arrival
            rayptmr2(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            
        end
        
        % azimuth from station to all tremors
        [~,az] = distance(stlasel2(i),stlosel2(i),tmrsel(:,7),tmrsel(:,6));
        aztmr2(:,i) = az;
    end
    toc
    save(strcat(datapath,'/ttsraypazi_tmr_bnd2_',prefix,'_',velflag,'.mat'),'ttstmr2','rayptmr2',...
         'aztmr2');
    
    % calculate the ray parameter between all stations and events
    ttsevt2 = zeros(nevt2, nsta2);
    raypevt2 = zeros(nevt2, nsta2);
    azevt2 = zeros(nevt2, nsta2);
    tic
    for i = 1: nsta2
        for j = 1: nevt2
            % travel time and ray p
            tt=tauptime('mod',velmod,'dep',evtbnd2(j,10),'ph','s,S','evt',...
                        [evtbnd2(j,9) evtbnd2(j,8)],'sta',[stlasel2(i) stlosel2(i)]);
            ttsevt2(j,i) = tt(1).time;   % choose the first S arrival
            raypevt2(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            
        end
        
        % azimuth from station to all tremors
        [~,az] = distance(stlasel2(i),stlosel2(i),evtbnd2(:,9),evtbnd2(:,8));
        azevt2(:,i) = az;
        
    end
    toc
    save(strcat(datapath,'/ttsraypazi_evt_bnd2_',prefix,'_',velflag,'.mat'),'ttsevt2','raypevt2',...
         'azevt2');
    
else
    data = load(strcat(datapath,'/ttsraypazi_tmr_bnd2_',prefix,'_',velflag,'.mat'));
    ttstmr2 = data.ttstmr2;
    rayptmr2 = data.rayptmr2;
    aztmr2 = data.aztmr2;
    
    data = load(strcat(datapath,'/ttsraypazi_evt_bnd2_',prefix,'_',velflag,'.mat'));
    ttsevt2 = data.ttsevt2;
    raypevt2 = data.raypevt2;
    azevt2 = data.azevt2;

end
%%%%%


%% visualization of ray parameter and azimuth distribution
%%%%% for bound 1 events and bound 0 tremors
f.fig=figure;
widin = 15;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow=6;
ncol=8;
for i = 1: nsta1
    f.ax(i) = subplot(nrow,ncol,i);
    hold(f.ax(i),'on');
    f.ax(i).Box='on';
    scatter(f.ax(i),rayptmr1(:,i),aztmr1(:,i),2,[0.3 0.3 0.3],'filled','o');
%     scatter(f.ax(i),raypevt1(:,i),azevt1(:,i),2,[0.5 0.8 1],'filled','o');
    text(f.ax(i),0.25,0.9,stnmsel1(i,:),'unit','normalized');
    hold(f.ax(i),'off');
        
end
xlabel(f.ax(nsta1),'Ray parameter (s/deg)');
ylabel(f.ax(nsta1),'Azimuth (deg)');
supertit(f.ax,strcat(prefix,{' bound 1'}));

% save figure
print(f.fig,'-depsc2',strcat(figpath,'/evt_tmr_rayp_azi_distri.',prefix,'.bd1.dep',...
      num2str(depran),'.eps'));
  
f.fig=figure;
widin = 15;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow=6;
ncol=8;
for i = 1: nsta1
    f.ax(i) = subplot(nrow,ncol,i);
    hold(f.ax(i),'on');
    f.ax(i).Box='on';
%     scatter(f.ax(i),rayptmr1(:,i),aztmr1(:,i),2,[0.3 0.3 0.3],'filled','o');
    scatter(f.ax(i),raypevt1(:,i),azevt1(:,i),2,[0.5 0.8 1],'filled','o');
    text(f.ax(i),0.25,0.9,stnmsel1(i,:),'unit','normalized');
    hold(f.ax(i),'off');
        
end
xlabel(f.ax(nsta1),'Ray parameter (s/deg)');
ylabel(f.ax(nsta1),'Azimuth (deg)');
supertit(f.ax,strcat(prefix,{' bound 1'}));

% save figure
print(f.fig,'-depsc2',strcat(figpath,'/evt_tmr_rayp_azi_distri.',prefix,'.bd1.dep',...
      num2str(depran),'.eps'));  
%%%%%

%%%%% for bound 2 events and bound 0 tremors
f2.fig=figure;
widin = 19;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow=6;
ncol=10;
for i = 1: nsta2
    f2.ax(i) = subplot(nrow,ncol,i);
    hold(f2.ax(i),'on');
    f2.ax(i).Box='on';
    scatter(f2.ax(i),rayptmr2(:,i),aztmr2(:,i),2,[0.3 0.3 0.3],'filled','o');
%     scatter(f2.ax(i),raypevt2(:,i),azevt2(:,i),2,[0.5 0.8 1],'filled','o');
    text(f2.ax(i),0.25,0.9,stnmsel2(i,:),'unit','normalized');
    hold(f2.ax(i),'off');
        
end
xlabel(f2.ax(nsta2),'Ray parameter (s/deg)');
ylabel(f2.ax(nsta2),'Azimuth (deg)');
supertit(f2.ax,strcat(prefix,{' bound 2'}));

% save figure
print(f2.fig,'-depsc2',strcat(figpath,'/evt_tmr_rayp_azi_distri.',prefix,'.bd2.dep',...
      num2str(depran),'.eps'));
  
f2.fig=figure;
widin = 19;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow=6;
ncol=10;
for i = 1: nsta2
    f2.ax(i) = subplot(nrow,ncol,i);
    hold(f2.ax(i),'on');
    f2.ax(i).Box='on';
%     scatter(f2.ax(i),rayptmr2(:,i),aztmr2(:,i),2,[0.3 0.3 0.3],'filled','o');
    scatter(f2.ax(i),raypevt2(:,i),azevt2(:,i),2,[0.5 0.8 1],'filled','o');
    text(f2.ax(i),0.25,0.9,stnmsel2(i,:),'unit','normalized');
    hold(f2.ax(i),'off');
        
end
xlabel(f2.ax(nsta2),'Ray parameter (s/deg)');
ylabel(f2.ax(nsta2),'Azimuth (deg)');
supertit(f2.ax,strcat(prefix,{' bound 2'}));

% save figure
print(f2.fig,'-depsc2',strcat(figpath,'/evt_tmr_rayp_azi_distri.',prefix,'.bd2.dep',...
      num2str(depran),'.eps'));  
%%%%%
  

%% set a tolerence to the difference in azimuth and ray p between tremors and events
% distribution of events should have no effects on the tolerence of azimuth and ray p
% in other words, the two tolerences are independent of events

%%%%%%%% for bound 1 events
%%% max azimuth difference, [deg]
aztol1 = 1*ones(nsta1,1);
% kfact = 2;
% aztol1 = kfact*((max(aztmr1)-min(aztmr1))./size(aztmr1,1));

%%% max ray p difference, [s/deg]
rayptol1 = 0.1*ones(nsta1,1);
% rayptol1 = kfact*((max(rayptmr1)-min(rayptmr1))./size(rayptmr1,1));
            
%%% relative time in days for the first day of tremor catalog           
d1 = datetime(tmrsel(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

recalflag = 0;

if recalflag
    tmrava1 = [];    % availble tremors
    evtava1 = [];    % paired availble events, one tremor could correspond to >=1 events
    tmravat1 = [];    % availble tremors, must at the same time, ie, within the same hour
    evtavat1 = [];    % paired availble events, one tremor could correspond to >=1 events
    tmravant1 = [];    % availble tremors
    evtavant1 = [];    % paired availble events, one tremor could correspond to >=1 events
    for i = 1: nsta1
        for j = 1: ntmr1  %
            % find all candidate events that has the 'same' ray p and azimuth as the tremor and pair
            % them, give them a label to keep tracking
            ind1 = find(abs(raypevt1(:,i)-rayptmr1(j,i))<=rayptol1(i) & ...
                abs(azevt1(:,i)-aztmr1(j,i))<=aztol1(i));
            if ~isempty(ind1)
                % add a label of tremor and station for tracking
                % so that the last 2 cols are flag to which tremor, and flag to which station
                tmrtmp = [tmrsel(j,:) j i];
                
                % add a label of tremor for tracking to available events as well
                evttmp = [evtbnd1(ind1,:) j*ones(length(ind1),1) i*ones(length(ind1),1)];
                
                % all available tremors and events
                tmrava1 = [tmrava1; tmrtmp];
                evtava1 = [evtava1; evttmp];
                
                % relative origin time in days of all tremors
                di = datetime(tmrtmp(1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
                ottmr = caldays(between(d1,di,'days')) + tmrtmp(5)/24;      % unit is day
                
                % relative origin time in days of all events
                di = datetime(evttmp(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
                otevt = caldays(between(d1,di,'days')) + (evttmp(:,5)+evttmp(:,6)/60+evttmp(:,7)...
                                /3600)/24;
                
                % travel time difference in day, this difference is theoretically in secs, so is
                % negligible compared with the tremor detection time resolution (hr)
                % So, this ttdiff is fine to delete
                ttdiff = (ttstmr1(j,i)-ttsevt1(ind1,i))./3600./24;    % must convert sec to day
                
                % find all pairs that must pass through the same area at the same time
                ind2 = find(otevt>=ottmr+ttdiff & otevt<ottmr+ttdiff+1/24);   % must within the same hr
                if ~isempty(ind2)
                    evttmp2 = evttmp(ind2, :);
                    tmravat1 = [tmravat1; tmrtmp];
                    evtavat1 = [evtavat1; evttmp2];
                    
                    % the remaining part is the ones pass through the same area NOT at the same time
                    tmravant1 = [tmravant1; tmrtmp];
                    tmp = evttmp;
                    tmp(ind2, :) = [];
                    evtavant1 = [evtavant1; tmp];
                else
                    tmravant1 = [tmravant1; tmrtmp];
                    evtavant1 = [evtavant1; evttmp];
                end
            end
        end
    end
    save(strcat(datapath,'/sameraypath_bnd1_',prefix,'_',velflag,'.mat'),'tmrava1','evtava1',...
         'tmravat1','evtavat1','tmravant1','evtavant1');
else
    data = load(strcat(datapath,'/sameraypath_bnd1_',prefix,'_',velflag,'.mat'));
    tmrava1 = data.tmrava1;
    evtava1 = data.evtava1;
    tmravat1 = data.tmravat1;
    evtavat1 = data.evtavat1;
    tmravant1 = data.tmravant1;
    evtavant1 = data.evtavant1;
end
%%%%%%%%

%%
%%%%%%%% for bound 2 events
%%% max azimuth difference, [deg]
% aztol = 0.05;
kfact = 2;
aztol2 = kfact*((max(aztmr2)-min(aztmr2))./size(aztmr2,1));

%%% max ray p difference, [s/rad]
% rayptol = 1;   
rayptol2 = kfact*((max(rayptmr2)-min(rayptmr2))./size(rayptmr2,1));
            
%%% relative time in days for the first day of tremor catalog           
d1 = datetime(tmrsel(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

recalflag = 0;

if recalflag
    tmrava2 = [];    % availble tremors
    evtava2 = [];    % paired availble events, one tremor could correspond to >=1 events
    tmravat2 = [];    % availble tremors, must at the same time, ie, within the same hour
    evtavat2 = [];    % paired availble events, one tremor could correspond to >=1 events
    tmravant2 = [];    % availble tremors
    evtavant2 = [];    % paired availble events, one tremor could correspond to >=1 events
    for i = 1: nsta2
        for j = 1: ntmr2  %
            % find all candidate events that has the 'same' ray p and azimuth as the tremor and pair
            % them, give them a label to keep tracking
            ind1 = find(abs(raypevt2(:,i)-rayptmr2(j,i))<=rayptol2(i) & ...
                abs(azevt2(:,i)-aztmr2(j,i))<=aztol2(i));
            if ~isempty(ind1)
                % add a label of tremor and station for tracking
                % so that the last 2 cols are flag to which tremor, and flag to which station
                tmrtmp = [tmrsel(j,:) j i];
                
                % add a label of tremor for tracking to available events as well
                evttmp = [evtbnd2(ind1,:) j*ones(length(ind1),1) i*ones(length(ind1),1)];
                
                % all available tremors and events
                tmrava2 = [tmrava2; tmrtmp];
                evtava2 = [evtava2; evttmp];
                
                % relative origin time in days of all tremors
                di = datetime(tmrtmp(1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
                ottmr = caldays(between(d1,di,'days')) + tmrtmp(5)/24;      % unit is day
                
                % relative origin time in days of all events
                di = datetime(evttmp(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
                otevt = caldays(between(d1,di,'days')) + (evttmp(:,5)+evttmp(:,6)/60+evttmp(:,6)/3600)/24;
                
                % travel time difference in day, this difference is theoretically in secs, so is
                % negligible compared with the tremor detection time resolution (hr)
                % So, this ttdiff is fine to delete
                ttdiff = (ttstmr2(j,i)-ttsevt2(ind1,i))./3600./24;    % must convert sec to day
                
                % find all pairs that must pass through the same area at the same time
                ind2 = find(otevt>=ottmr+ttdiff & otevt<ottmr+ttdiff+1/24);   % must within the same hr
                if ~isempty(ind2)
                    evttmp2 = evttmp(ind2, :);
                    tmravat2 = [tmravat2; tmrtmp];
                    evtavat2 = [evtavat2; evttmp2];
                    
                    % the remaining part is the ones pass through the same area NOT at the same time
                    tmravant2 = [tmravant2; tmrtmp];
                    tmp = evttmp;
                    tmp(ind2, :) = [];
                    evtavant2 = [evtavant2; tmp];
                else
                    tmravant2 = [tmravant2; tmrtmp];
                    evtavant2 = [evtavant2; evttmp];
                end
            end
        end
    end
    save(strcat(datapath,'/sameraypath_bnd2_',prefix,'.mat'),'tmrava2','evtava2',...
         'tmravat2','evtavat2','tmravant2','evtavant2');
else
    data = load(strcat(datapath,'/sameraypath_bnd2_',prefix,'.mat'));
    tmrava2 = data.tmrava2;
    evtava2 = data.evtava2;
    tmravat2 = data.tmravat2;
    evtavat2 = data.evtavat2;
    tmravant2 = data.tmravant2;
    evtavant2 = data.evtavant2;
end    
%%%%%%%%


%% visualization of spatial distribtuion of those share the same path and occur in the same time
f3.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f3.ax = gca;
hold(f3.ax,'on');
grid(f3.ax, 'on');
f3.ax.Box='on';
% if regflag == 1  % means western shikoku
%     plot(f3.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
%          [0.6 0.6 0.6],'linew',2);
% elseif regflag == 2  % means kii pen
%     plot(f3.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
%         [0.6 0.6 0.6],'linew',2);
% end
axis(f3.ax, 'equal');
if regflag == 1
    xlim(f3.ax,[131.5 134.5]);
    ylim(f3.ax,[33 35]);
elseif regflag == 2
    xlim(f3.ax,[135 136.7]);
    ylim(f3.ax,[33.5 35]);
end

coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(f3.ax,coast(:,1),coast(:,2),'black','linew',0.5);

lat = 31:0.01:36;
lon = 131:0.01:138;
[longrd, latgrd] = meshgrid(lon,lat);
depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');
contour(f3.ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.8 0.8 0.8]);

%%%%%% for bound 1 events and tremors
indsta = unique(tmravat1(:,end));
for i = 1: length(indsta)
    scatter(f3.ax,stlosel1(indsta(i)), stlasel1(indsta(i)),70,'y','filled','^','markeredgec','k');
    text(f3.ax,stlosel1(indsta(i)), stlasel1(indsta(i))+0.03, stnmsel1(indsta(i),:), 'fontsize',8);
    tmrtmp = tmravat1(tmravat1(:,end)==indsta(i),:);
    p1=scatter(f3.ax,tmrtmp(:,6),tmrtmp(:,7),30,'r','filled','o','markeredgec','k');
    evttmp = evtavat1(evtavat1(:,end)==indsta(i),:);
    p2=scatter(f3.ax,evttmp(:,8),evttmp(:,9),30,'b','filled','o','markeredgec','k');
    for j= 1: size(tmrtmp)
        plot(f3.ax,[stlosel1(indsta(i)) tmrtmp(j,6)], [stlasel1(indsta(i)) tmrtmp(j,7)],'color',...
             [0.6 0.6 0.6]);
    end
    for j= 1: size(evttmp)
        plot(f3.ax,[stlosel1(indsta(i)) evttmp(j,8)], [stlasel1(indsta(i)) evttmp(j,9)],'color',...
             [0.6 0.6 0.6]);
    end
end
%%%%%%

% %%%%%% for bound 2 events and tremors
% indsta = unique(tmravat2(:,end));
% for i = 1: length(indsta)
%     scatter(f3.ax,stlosel2(indsta(i)), stlasel2(indsta(i)),70,'y','filled','^','markeredgec','k');
%     text(f3.ax,stlosel2(indsta(i)), stlasel2(indsta(i))+0.03, stnmsel2(indsta(i),:), 'fontsize',8);
%     tmrtmp = tmravat2(tmravat2(:,end)==indsta(i),:);
%     scatter(f3.ax,tmrtmp(:,6),tmrtmp(:,7),30,'r','filled','o','markeredgec','k');
%     evttmp = evtavat2(evtavat2(:,end)==indsta(i),:);
%     scatter(f3.ax,evttmp(:,8),evttmp(:,9),30,'b','filled','o','markeredgec','k');
%     for j= 1: size(tmrtmp)
%         plot(f3.ax,[stlosel2(indsta(i)) tmrtmp(j,6)], [stlasel2(indsta(i)) tmrtmp(j,7)],'color',...
%              [0.6 0.6 0.6]);
%     end
%     for j= 1: size(evttmp)
%         plot(f3.ax,[stlosel2(indsta(i)) evttmp(j,8)], [stlasel2(indsta(i)) evttmp(j,9)],'color',...
%              [0.6 0.6 0.6]);
%     end
% end
% %%%%%%

legend([p1,p2],{'Tremors','Events'},'location','northwest','fontsize',8);
title(f3.ax,prefix);
hold(f3.ax,'off');

% save figure
print(f3.fig,'-dpdf',strcat(figpath,'/evt_tmr_samerayp_sametime.',prefix,'.dep',...
      num2str(depran),'_',velflag,'.pdf'));
  
  
  
%% visualization of frequency-magnitude distribution that share the same ray path
f4.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f4.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow = 2;
ncol = 2;
for i = 1: ncol*nrow
    f4.ax(i) = subplot(nrow,ncol,i);
end

binw = 0.5;
f4.ax(1) = pltevtmagfreq(f4.ax(1),evtavat1(:,11),binw);
f4.ax(2) = pltevtmagfreq(f4.ax(2),evtavant1(:,11),binw);
f4.ax(3) = pltevtmagfreq(f4.ax(3),evtavat2(:,11),binw);
f4.ax(4) = pltevtmagfreq(f4.ax(4),evtavant2(:,11),binw);

hold(f4.ax(1),'on');
text(f4.ax(1),0.5,0.2,{'close in time'},'unit','normalized');
text(f4.ax(1),0.5,0.1,{'bound 1'},'unit','normalized');
xlim(f4.ax(1),[-1 5]);
hold(f4.ax(1),'off');

hold(f4.ax(2),'on');
text(f4.ax(2),0.5,0.2,{'not close in time'},'unit','normalized');
xlim(f4.ax(2),[-1 5]);
hold(f4.ax(2),'off');

hold(f4.ax(3),'on');
text(f4.ax(3),0.5,0.1,{'bound 2'},'unit','normalized');
xlim(f4.ax(3),[-1 5]);
hold(f4.ax(3),'off');

xlim(f4.ax(4),[-1 5]);

supertit(f4.ax,prefix);        

% save figure
print(f4.fig,'-dpdf',strcat(figpath,'/evt_freq_mag_samerayp.',prefix,'.dep',...
      num2str(depran),'.pdf'));
  
  
  
%% calculate the event rate during and outside tremor times

%%%%%% for bound 0 (the rectangle) events and tremors
datemin = min(obara(:,1));   % in bound 0, tremor and events have same time coverage
datemax = max(obara(:,1));

datetmr = unique(obara(:,1));
lentmr0 = length(datetmr);

d1 = datetime(datemin,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
d2 = datetime(datemax,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
lenall0 = caldays(between(d1,d2,'days'))+1;

ind = [];
for i = 1: lentmr0
   indtmp = find(datetmr(i) == evtall(:,1));
   ind = [ind; indtmp];
end
evtintmr0 = evtall(ind, :);
nintmr0 = size(evtintmr0,1);

indall = 1: size(evtall,1);
ind2 = setdiff(indall,ind);

evtouttmr0 = evtall(ind2, :);
nouttmr0 = size(evtouttmr0,1);

rtintmr0 = nintmr0/lentmr0;     % num per day
rtouttmr0 = nouttmr0/(lenall0-lentmr0);     
  

%%%%%% for bound 1 events and tremors
datemin = max(min(tmrbnd1(:,1)), min(evtbnd1(:,1)));   
datemax = min(max(tmrbnd1(:,1)), max(evtbnd1(:,1)));

tmrtmp = tmrbnd1(tmrbnd1(:,1)>=datemin & tmrbnd1(:,1)<=datemax, :);
evttmp = evtbnd1(evtbnd1(:,1)>=datemin & evtbnd1(:,1)<=datemax, :);

datetmr = unique(tmrtmp(:,1));
lentmr1 = length(datetmr);

d1 = datetime(datemin,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
d2 = datetime(datemax,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
lenall1 = caldays(between(d1,d2,'days'))+1;

ind = [];
for i = 1: lentmr1
   indtmp = find(datetmr(i) == evttmp(:,1));
   ind = [ind; indtmp];
end
evtintmr1 = evttmp(ind, :);
nintmr1 = size(evtintmr1,1);

indall = 1: size(evttmp,1);
ind2 = setdiff(indall,ind);

evtouttmr1 = evttmp(ind2, :);
nouttmr1 = size(evtouttmr1,1);

rtintmr1 = nintmr1/lentmr1;     % num per day
rtouttmr1 = nouttmr1/(lenall1-lentmr1);


%%%%%% for bound 2 events and tremors
datemin = max(min(tmrbnd2(:,1)), min(evtbnd2(:,1)));   
datemax = min(max(tmrbnd2(:,1)), max(evtbnd2(:,1)));

tmrtmp = tmrbnd2(tmrbnd2(:,1)>=datemin & tmrbnd2(:,1)<=datemax, :);
evttmp = evtbnd2(evtbnd2(:,1)>=datemin & evtbnd2(:,1)<=datemax, :);

datetmr = unique(tmrtmp(:,1));
lentmr2 = length(datetmr);

d1 = datetime(datemin,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
d2 = datetime(datemax,'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
lenall2 = caldays(between(d1,d2,'days'))+1;

ind = [];
for i = 1: lentmr2
   indtmp = find(datetmr(i) == evttmp(:,1));
   ind = [ind; indtmp];
end
evtintmr2 = evttmp(ind, :);
nintmr2 = size(evtintmr2,1);

indall = 1: size(evttmp,1);
ind2 = setdiff(indall,ind);

evtouttmr2 = evttmp(ind2, :);
nouttmr2 = size(evtouttmr2,1);

rtintmr2 = nintmr2/lentmr2;
rtouttmr2 = nouttmr2/(lenall2-lentmr2);
        
        
% plot the mag-freq chart and see if there is any difference         
f5.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow = 3;
ncol = 2;
for i = 1: ncol*nrow
    f5.ax(i) = subplot(nrow,ncol,i);
end

binw = 0.5;
[f5.ax(1),h1] = pltevtmagfreq(f5.ax(1),evtintmr0(:,11),binw);
[f5.ax(2),h2] = pltevtmagfreq(f5.ax(2),evtouttmr0(:,11),binw);
[f5.ax(3),h3] = pltevtmagfreq(f5.ax(3),evtintmr1(:,11),binw);
[f5.ax(4),h4] = pltevtmagfreq(f5.ax(4),evtouttmr1(:,11),binw);
[f5.ax(5),h5] = pltevtmagfreq(f5.ax(5),evtintmr2(:,11),binw);
[f5.ax(6),h6] = pltevtmagfreq(f5.ax(6),evtouttmr2(:,11),binw);

rtintmr0m = h1.Values/lentmr0;
rtouttmr0m = h2.Values/(lenall0-lentmr0);
rtintmr1m = h3.Values/lentmr1;
rtouttmr1m = h4.Values/(lenall1-lentmr1);
rtintmr2m = h5.Values/lentmr2;
rtouttmr2m = h6.Values/(lenall2-lentmr2);

hold(f5.ax(1),'on');
text(f5.ax(1),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtintmr0)),'unit',...
     'normalized');
for i = 1: length(h1.Values)
    text(f5.ax(1),h1.BinEdges(i),h1.Values(i)+f5.ax(1).YLim(2)/30,sprintf('%.1f',rtintmr0m(i)),...
         'fontsize',8); 
end
text(f5.ax(1),0.65,0.3,{'in bound 0'},'unit','normalized');
xlim(f5.ax(1),[-1 6]);
title(f5.ax(1),'During tremor days');
hold(f5.ax(1),'off');        
        
hold(f5.ax(2),'on');
text(f5.ax(2),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtouttmr0)),'unit',...
     'normalized');
for i = 1: length(h2.Values)
    text(f5.ax(2),h2.BinEdges(i),h2.Values(i)+f5.ax(2).YLim(2)/30,sprintf('%.1f',rtouttmr0m(i)),...
         'fontsize',8); 
end
xlim(f5.ax(2),[-1 6]);
title(f5.ax(2),'Outside tremor days');
hold(f5.ax(2),'off'); 

hold(f5.ax(3),'on');
text(f5.ax(3),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtintmr1)),'unit',...
     'normalized');
for i = 1: length(h3.Values)
    text(f5.ax(3),h3.BinEdges(i),h3.Values(i)+f5.ax(3).YLim(2)/30,sprintf('%.1f',rtintmr1m(i)),...
         'fontsize',8); 
end
text(f5.ax(3),0.65,0.3,{'in bound 1'},'unit','normalized');
xlim(f5.ax(3),[-1 6]);
hold(f5.ax(3),'off'); 

hold(f5.ax(4),'on');
text(f5.ax(4),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtouttmr1)),'unit',...
     'normalized');
for i = 1: length(h4.Values)
    text(f5.ax(4),h4.BinEdges(i),h4.Values(i)+f5.ax(4).YLim(2)/30,sprintf('%.1f',rtouttmr1m(i)),...
         'fontsize',8); 
end
xlim(f5.ax(4),[-1 6]);
hold(f5.ax(4),'off'); 

hold(f5.ax(5),'on');
text(f5.ax(5),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtintmr2)),'unit',...
     'normalized');
for i = 1: length(h5.Values)
    text(f5.ax(5),h5.BinEdges(i),h5.Values(i)+f5.ax(5).YLim(2)/30,sprintf('%.1f',rtintmr2m(i)),...
         'fontsize',8); 
end
text(f5.ax(5),0.65,0.3,{'in bound 2'},'unit','normalized');
xlim(f5.ax(5),[-1 6]);
hold(f5.ax(5),'off'); 

hold(f5.ax(6),'on');
text(f5.ax(6),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtouttmr2)),'unit',...
     'normalized');
for i = 1: length(h6.Values)
    text(f5.ax(6),h6.BinEdges(i),h6.Values(i)+f5.ax(6).YLim(2)/30,sprintf('%.1f',rtouttmr2m(i)),...
         'fontsize',8); 
end
xlim(f5.ax(6),[-1 6]);
hold(f5.ax(6),'off');  

supertit(f5.ax,prefix);

% save figure
print(f5.fig,'-dpdf',strcat(figpath,'/evt_freq_mag_diff_bounds.',prefix,'.dep',...
      num2str(depran),'.pdf'));   




























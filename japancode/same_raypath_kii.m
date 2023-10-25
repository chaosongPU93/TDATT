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


% read the regular events catalog
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

%%% check the event catalog more carefully, retain the natural EQ only
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


%% set the further boundary to regular events 
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

% because a large region in the SE part of the rectangle is blank, thus redundant
%%% bound type 1, a smaller rectangle
if regflag == 1  % means western shikoku
    nrectx = [134.1 134.7 132.11 131.51]';
    nrecty = [35.1 33.7 32.59 33.99]';
elseif regflag == 2  % means kii pen
    nrectx = [134.6 135.2 136.22 135.62]';
    nrecty = [34.1 33.2 33.88 34.78]';
end
[ind,areaall] = boundary(nrectx,nrecty,0);
bnd1 = [nrectx(ind) nrecty(ind)];

%%% bound type 2, irregular shape hand picked by choosing the co-located region as tremors
bnd2 = [135.8353   34.2355;
        135.9033   34.1033;
        135.8125   33.9557;
        135.7605   33.8423;
        135.6840   33.7305;
        135.5730   33.6992;
        135.4980   33.7128;
        135.4265   33.7350;
        135.4012   33.7782;
        135.3295   33.7925;
        135.3128   33.8605;
        135.3430   33.9337;
        135.4308   33.9822;
        135.5017   34.0917;
        135.6038   34.2657;
        135.8353   34.2355];

% events & tremor inside bound 1
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),bnd1(:,1),bnd1(:,2));
isinbnd1 = is | ion;
evtbnd1 = evtall(isinbnd1 == 1, :);
[is,ion] = inpolygon(obara(:,6),obara(:,7),bnd1(:,1),bnd1(:,2));
isinbnd1 = is | ion;
tmrbnd1 = obara(isinbnd1 == 1, :);

% events & tremor inside bound 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),bnd2(:,1),bnd2(:,2));
isinbnd2 = is | ion;
evtbnd2 = evtall(isinbnd2 == 1, :);
[is,ion] = inpolygon(obara(:,6),obara(:,7),bnd2(:,1),bnd2(:,2));
isinbnd2 = is | ion;
tmrbnd2 = obara(isinbnd2 == 1, :);


%% create a cross-scetion profile, project all events & tremor onto it
% csp1x = linspace(135.8,135.2,1000);
% csp1y = linspace(33.6,34.5,1000);
% csp1 = [csp1x' csp1y'];
% 
% csp2x = linspace(135.53,135.4,1000);
% csp2y = linspace(33.42,34.63,1000);
% csp2 = [csp2x' csp2y'];

% cross-section region (rect)
csp = [135.5  33.4;
       135.38 34.6;
       135.2904 34.5910;
       135.4104 33.3910;
       135.5  33.4];
   
% the middle line of the cross-section
cspmidx = linspace(0.5*(csp(4,1)+csp(5,1)),0.5*(csp(2,1)+csp(3,1)),1000);
cspmidy = linspace(0.5*(csp(4,2)+csp(5,2)),0.5*(csp(2,2)+csp(3,2)),1000);
cspmid = [cspmidx' cspmidy'];

% assume tremors are on the slab interface, interpolate to get their depth
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
obara(:,8) = F(obara(:,6),obara(:,7));

% define the origin as the mid lower right point
lon0 = cspmid(1,1);
lat0 = cspmid(1,2);

% convert the cross-section profile to relative distance 
[cspmidr(:,1),cspmidr(:,2)] = absloc2relaloc(cspmid(:,1),cspmid(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(cspmidr(end,2),cspmidr(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
cspmidr(:,3) = F(cspmid(:,1),cspmid(:,2));
% rotate the coordinates to profile axis 
[cspmidr(:,4),cspmidr(:,5)] = coordinate_rot(cspmidr(:,1),cspmidr(:,2),rotang,0,0);

% convert events' locations to relative distance
[dx,dy] = absloc2relaloc(evtbnd1(:,8),evtbnd1(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtbnd1(:,13),evtbnd1(:,14)] = coordinate_rot(dx,dy,rotang,0,0);
% do the same to bound type 2 events
[dx,dy] = absloc2relaloc(evtbnd2(:,8),evtbnd2(:,9),lon0,lat0);
[evtbnd2(:,13),evtbnd2(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(obara(:,6),obara(:,7),lon0,lat0);
[obara(:,9),obara(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stxr,styr] = coordinate_rot(dx,dy,rotang,0,0);


%% keep only staions that can sample the portion above the tremor, rather than updip slab
%   Find the upper depth limit of tremors, any stations whose along profile distance
%	smaller than that at this depth limit can be discarded, because those stations only sample
%	the updip portion of the slab, which is definitely not the zone we are seeking.
[mindep, ind1] = min(obara(:,8));
[mindist, ind2] = min(obara(:,10));   % the two indexes might be the same, using either is fine. 

[is,ion] = inpolygon(stlo,stla,bnd0(:,1),bnd0(:,2));
isinbnd0 = is | ion;

indsel = find(styr>=mindist & isinbnd0==1);
stnmsel = stnm(indsel,:);    % station name, string
stlosel = stlo(indsel);    % station longitude
stlasel = stla(indsel);    % station latitude
stelsel = stel(indsel);    % station elevation

%% massive calculation
% time for travel time and ray p for all stations and tremors in bound 0 take 1610 s, ~0.5 h
% time for travel time and ray p for all stations and events in bound 1 take 3574 s, ~1 h
evtsel = evtbnd1;
tmrsel = obara;
nsta = size(stnmsel,1);
ntmr = size(tmrsel,1);
nevt = size(evtsel,1);

% recalculation flag
recalflag = 0;

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

if recalflag
    % calculate the ray parameter between all stations and tremors
    ttstmr = zeros(ntmr, nsta);
    rayptmr = zeros(ntmr, nsta);
    aztmr = zeros(ntmr, nsta);
    tic
    for i = 1: nsta
        for j = 1: ntmr
            % travel time and ray p
            tt=tauptime('mod',velmod,'dep',tmrsel(j,8),'ph','s,S','evt',[tmrsel(j,7) tmrsel(j,6)],...
                        'sta',[stlasel(i) stlosel(i)]);
            ttstmr(j,i) = tt(1).time;   % choose the first S arrival
            rayptmr(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            
        end
        
        % azimuth from station to all tremors
        [~,az] = distance(stlasel(i),stlosel(i),tmrsel(:,7),tmrsel(:,6));
        aztmr(:,i) = az;    % unit in [deg]
    end
    toc
    save(strcat(datapath,'/ttsraypazi_tmr_bnd0_',prefix,'_',velflag,'.mat'),'ttstmr','rayptmr',...
         'aztmr');
    
    % calculate the ray parameter between all stations and events
    ttsevt = zeros(nevt, nsta);
    raypevt = zeros(nevt, nsta);
    azevt = zeros(nevt, nsta);
    tic
    for i = 1: nsta
        for j = 1: nevt
            % travel time and ray p
            tt=tauptime('mod',velmod,'dep',evtsel(j,10),'ph','s,S','evt',[evtsel(j,9) evtsel(j,8)],...
                        'sta',[stlasel(i) stlosel(i)]);
            ttsevt(j,i) = tt(1).time;   % choose the first S arrival
            raypevt(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            
        end
        
        % azimuth from station to all tremors
        [~,az] = distance(stlasel(i),stlosel(i),evtsel(:,9),evtsel(:,8));
        azevt(:,i) = az;    % unit in [deg]
        
    end
    toc
    save(strcat(datapath,'/ttsraypazi_evt_bnd1_',prefix,'_',velflag,'.mat'),'ttsevt','raypevt',...
         'azevt');
    
else
    data = load(strcat(datapath,'/ttsraypazi_tmr_bnd0_',prefix,'_',velflag,'.mat'));
    ttstmr = data.ttstmr;
    rayptmr = data.rayptmr;
    aztmr = data.aztmr;
    
    data = load(strcat(datapath,'/ttsraypazi_evt_bnd1_',prefix,'_',velflag,'.mat'));
    ttsevt = data.ttsevt;
    raypevt = data.raypevt;
    azevt = data.azevt;

end


%% visualization of ray parameter and azimuth distribution
% how to make distribution of events and tremors stand out separately, since there are also overlaps
% between the two
%   1. plot one of them separately on two figures
%   2. plot them on top of each other in turn on two figures
%   3. plot still in one figure, use transparency, but not all stations on the same figure, will be
%      too large to recognize the transparency

f.fig=figure;
widin = 10;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow=5;
ncol=5;
for i = 1: nsta
    f.ax(i) = subplot(nrow,ncol,i);
    hold(f.ax(i),'on');
    f.ax(i).Box='on';
    grid(f.ax(i),'on');
    scatter(f.ax(i),rayptmr(:,i),aztmr(:,i),2,[0.3 0.3 0.3],'filled','o');
%     scatter(f.ax(i),raypevt(:,i),azevt(:,i),2,[0.5 0.8 1],'filled','o');
    text(f.ax(i),0.25,0.9,stnmsel(i,:),'unit','normalized');
    hold(f.ax(i),'off');
        
end
xlabel(f.ax(nsta),'Ray parameter (s/rad)');
ylabel(f.ax(nsta),'Azimuth (deg)');
supertit(f.ax,prefix);

% save figure
print(f.fig,'-depsc2',strcat(figpath,'/evt_tmr_rayp_azi_distri.',prefix,'.dep',...
      num2str(depran),'.eps'));

  
f.fig=figure;
widin = 10;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow=5;
ncol=5;
for i = 1: nsta
    f.ax(i) = subplot(nrow,ncol,i);
    hold(f.ax(i),'on');
    f.ax(i).Box='on';
    grid(f.ax(i),'on');
%     scatter(f.ax(i),rayptmr(:,i),aztmr(:,i),2,[0.3 0.3 0.3],'filled','o');
    scatter(f.ax(i),raypevt(:,i),azevt(:,i),2,[0.5 0.8 1],'filled','o');
    text(f.ax(i),0.25,0.9,stnmsel(i,:),'unit','normalized');
    hold(f.ax(i),'off');
        
end
xlabel(f.ax(nsta),'Ray parameter (s/rad)');
ylabel(f.ax(nsta),'Azimuth (deg)');
supertit(f.ax,prefix);

% save figure
print(f.fig,'-depsc2',strcat(figpath,'/evt_tmr_rayp_azi_distri.',prefix,'.dep',...
      num2str(depran),'.eps'));  
        
%% set a tolerence to the difference in azimuth and ray p between tremors and events
% distribution of events should have no effects on the tolerence of azimuth and ray p
% in other words, the two tolerences are independent of events

%%% max azimuth difference, [deg]
aztol = 1*ones(nsta,1);
% kfact = 2;
% aztol = kfact*((max(aztmr)-min(aztmr))./size(aztmr,1));

%%% max ray p difference, [s/deg]
rayptol = 0.1*ones(nsta,1);   
% rayptol = kfact*((max(rayptmr)-min(rayptmr))./size(rayptmr,1));
            
%%% relative time in days for the first day of tremor catalog           
d1 = datetime(tmrsel(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

recalflag = 1;

if recalflag
    tmrava = [];    % availble tremors
    evtava = [];    % paired availble events, one tremor could correspond to >=1 events
    tmravat = [];    % availble tremors, must at the same time, ie, within the same hour
    evtavat = [];    % paired availble events, one tremor could correspond to >=1 events
    tmravant = [];    % availble tremors
    evtavant = [];    % paired availble events, one tremor could correspond to >=1 events
    
    for i = 1: nsta
        for j = 1: ntmr  %
            % find all candidate events that has the 'same' ray p and azimuth as the tremor and pair
            % them, give them a label to keep tracking
            ind1 = find(abs(raypevt(:,i)-rayptmr(j,i))<=rayptol(i) & ...
                abs(azevt(:,i)-aztmr(j,i))<=aztol(i));
            if ~isempty(ind1)
                % add a label of tremor and station for tracking
                % so that the last 2 cols are flag to which tremor, and flag to which station
                tmrtmp = [tmrsel(j,:) j i];
                
                % add a label of tremor for tracking to available events as well
                evttmp = [evtsel(ind1,:) j*ones(length(ind1),1) i*ones(length(ind1),1)];
                
                % all available tremors and events
                tmrava = [tmrava; tmrtmp];
                evtava = [evtava; evttmp];
                
                % relative origin time in days of all tremors
                di = datetime(tmrtmp(1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
                ottmr = caldays(between(d1,di,'days')) + tmrtmp(5)/24;      % unit is day
                
                % relative origin time in days of all events
                di = datetime(evttmp(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
                otevt = caldays(between(d1,di,'days')) + (evttmp(:,5)+evttmp(:,6)/60+evttmp(:,7)/3600)/24;
                
                % travel time difference in day, this difference is theoretically in secs, so is
                % negligible compared with the tremor detection time resolution (hr)
                % So, this ttdiff is fine to delete
                ttdiff = (ttstmr(j,i)-ttsevt(ind1,i))./3600./24;    % must convert sec to day
                
                % find all pairs that must pass through the same area at the same time
                ind2 = find(otevt>=ottmr+ttdiff & otevt<ottmr+ttdiff+1/24);   % must within the same hr
                if ~isempty(ind2)
                    evttmp2 = evttmp(ind2, :);
                    tmravat = [tmravat; tmrtmp];
                    evtavat = [evtavat; evttmp2];
                    
                    % the remaining part is the ones pass through the same area NOT at the same time
                    tmravant = [tmravant; tmrtmp];
                    tmp = evttmp;
                    tmp(ind2, :) = [];
                    evtavant = [evtavant; tmp];
                else
                    tmravant = [tmravant; tmrtmp];
                    evtavant = [evtavant; evttmp];
                end
            end
        end
    end
    save(strcat(datapath,'/sameraypath_bnd1_',prefix,'.mat'),'tmrava','evtava','tmravat',...
         'evtavat','tmravant','evtavant');
else
    data = load(strcat(datapath,'/sameraypath_bnd1_',prefix,'.mat'));
    tmrava = data.tmrava;
    evtava = data.evtava;
    tmravat = data.tmravat;
    evtavat = data.evtavat;
    tmravant = data.tmravant;
    evtavant = data.evtavant;
end


%% visualization of spatial distribtuion of those share the same path and occur in the same time
f2.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f2.ax = gca;
hold(f2.ax,'on');
grid(f2.ax, 'on');
f2.ax.Box='on';
% if regflag == 1  % means western shikoku
%     plot(f2.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
%          [0.6 0.6 0.6],'linew',2);
% elseif regflag == 2  % means kii pen
%     plot(f2.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
%         [0.6 0.6 0.6],'linew',2);
% end
axis(f2.ax, 'equal');
if regflag == 1
    xlim(f2.ax,[131 135]);
    ylim(f2.ax,[32.3 35.2]);
elseif regflag == 2
    xlim(f2.ax,[135 136.7]);
    ylim(f2.ax,[33.5 35]);
end

coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(f2.ax,coast(:,1),coast(:,2),'black','linew',0.5);

lat = 31:0.01:36;
lon = 131:0.01:138;
[longrd, latgrd] = meshgrid(lon,lat);
depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');
contour(f2.ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.8 0.8 0.8]);

indsta = unique(tmravat(:,end));
for i = 1: length(indsta)
    scatter(f2.ax,stlosel(indsta(i)), stlasel(indsta(i)),70,'y','filled','^','markeredgec','k');
    text(f2.ax,stlosel(indsta(i)), stlasel(indsta(i))+0.03, stnmsel(indsta(i),:), 'fontsize',8);
    tmrtmp = tmravat(tmravat(:,end)==indsta(i),:);
    p1=scatter(f2.ax,tmrtmp(:,6),tmrtmp(:,7),30,'r','filled','o','markeredgec','k');
    evttmp = evtavat(evtavat(:,end)==indsta(i),:);
    p2=scatter(f2.ax,evttmp(:,8),evttmp(:,9),30,'b','filled','o','markeredgec','k');
    for j= 1: size(tmrtmp)
        plot(f2.ax,[stlosel(indsta(i)) tmrtmp(j,6)], [stlasel(indsta(i)) tmrtmp(j,7)],'color',...
             [0.6 0.6 0.6]);
    end
    for j= 1: size(evttmp)
        plot(f2.ax,[stlosel(indsta(i)) evttmp(j,8)], [stlasel(indsta(i)) evttmp(j,9)],'color',...
             [0.6 0.6 0.6]);
    end
end
legend([p1,p2],{'Tremors','Events'},'location','northwest','fontsize',8);
title(f2.ax,prefix);
hold(f2.ax,'off');

% save figure
print(f2.fig,'-dpdf',strcat(figpath,'/evt_tmr_samerayp_sametime.',prefix,'.dep',...
      num2str(depran),'.pdf'));
  

%% visualization of frequency-magnitude distribution that share the same ray path
f3.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow = 1;
ncol = 2;
for i = 1: ncol*nrow
    f3.ax(i) = subplot(nrow,ncol,i);
end

binw = 0.5;
f3.ax(1) = pltevtmagfreq(f3.ax(1),evtavat(:,11),binw);
f3.ax(2) = pltevtmagfreq(f3.ax(2),evtavant(:,11),binw);

hold(f3.ax(1),'on');
text(f3.ax(1),0.5,0.2,{'close in time'},'unit','normalized');
xlim(f3.ax(1),[-1 5]);
hold(f3.ax(1),'off');

hold(f3.ax(2),'on');
text(f3.ax(2),0.5,0.2,{'not close in time'},'unit','normalized');
xlim(f3.ax(2),[-1 5]);
hold(f3.ax(2),'off');

supertit(f3.ax,prefix);        

% save figure
print(f3.fig,'-dpdf',strcat(figpath,'/evt_freq_mag_samerayp.',prefix,'.dep',...
      num2str(depran),'.pdf'));

 
  
%% calculate the event rate during and outside tremor times
%%% 1. events and tremors in bound 0
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

rtintmr0 = nintmr0/lentmr0;
rtouttmr0 = nouttmr0/(lenall0-lentmr0);
  

%%% 2. events and tremors in bound 1
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

rtintmr1 = nintmr1/lentmr1;
rtouttmr1 = nouttmr1/(lenall1-lentmr1);


%%% 3. events and tremors in bound 2        
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
f4.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f4.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

nrow = 3;
ncol = 2;
for i = 1: ncol*nrow
    f4.ax(i) = subplot(nrow,ncol,i);
end

binw = 0.5;
[f4.ax(1),h1] = pltevtmagfreq(f4.ax(1),evtintmr0(:,11),binw);
[f4.ax(2),h2] = pltevtmagfreq(f4.ax(2),evtouttmr0(:,11),binw);
[f4.ax(3),h3] = pltevtmagfreq(f4.ax(3),evtintmr1(:,11),binw);
[f4.ax(4),h4] = pltevtmagfreq(f4.ax(4),evtouttmr1(:,11),binw);
[f4.ax(5),h5] = pltevtmagfreq(f4.ax(5),evtintmr2(:,11),binw);
[f4.ax(6),h6] = pltevtmagfreq(f4.ax(6),evtouttmr2(:,11),binw);

rtintmr0m = h1.Values/lentmr0;
rtouttmr0m = h2.Values/(lenall0-lentmr0);
rtintmr1m = h3.Values/lentmr1;
rtouttmr1m = h4.Values/(lenall1-lentmr1);
rtintmr2m = h5.Values/lentmr2;
rtouttmr2m = h6.Values/(lenall2-lentmr2);

hold(f4.ax(1),'on');
text(f4.ax(1),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtintmr0)),'unit',...
     'normalized');
for i = 1: length(h1.Values)
    text(f4.ax(1),h1.BinEdges(i),h1.Values(i)+f4.ax(1).YLim(2)/30,sprintf('%.1f',rtintmr0m(i)),...
         'fontsize',8); 
end
text(f4.ax(1),0.65,0.3,{'in bound 0'},'unit','normalized');
xlim(f4.ax(1),[-1 6]);
title(f4.ax(1),'During tremor days');
hold(f4.ax(1),'off');        
        
hold(f4.ax(2),'on');
text(f4.ax(2),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtouttmr0)),'unit',...
     'normalized');
for i = 1: length(h2.Values)
    text(f4.ax(2),h2.BinEdges(i),h2.Values(i)+f4.ax(2).YLim(2)/30,sprintf('%.1f',rtouttmr0m(i)),...
         'fontsize',8); 
end
xlim(f4.ax(2),[-1 6]);
title(f4.ax(2),'Outside tremor days');
hold(f4.ax(2),'off'); 

hold(f4.ax(3),'on');
text(f4.ax(3),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtintmr1)),'unit',...
     'normalized');
for i = 1: length(h3.Values)
    text(f4.ax(3),h3.BinEdges(i),h3.Values(i)+f4.ax(3).YLim(2)/30,sprintf('%.1f',rtintmr1m(i)),...
         'fontsize',8); 
end
text(f4.ax(3),0.65,0.3,{'in bound 1'},'unit','normalized');
xlim(f4.ax(3),[-1 6]);
hold(f4.ax(3),'off'); 

hold(f4.ax(4),'on');
text(f4.ax(4),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtouttmr1)),'unit',...
     'normalized');
for i = 1: length(h4.Values)
    text(f4.ax(4),h4.BinEdges(i),h4.Values(i)+f4.ax(4).YLim(2)/30,sprintf('%.1f',rtouttmr1m(i)),...
         'fontsize',8); 
end
xlim(f4.ax(4),[-1 6]);
hold(f4.ax(4),'off'); 

hold(f4.ax(5),'on');
text(f4.ax(5),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtintmr2)),'unit',...
     'normalized');
for i = 1: length(h5.Values)
    text(f4.ax(5),h5.BinEdges(i),h5.Values(i)+f4.ax(5).YLim(2)/30,sprintf('%.1f',rtintmr2m(i)),...
         'fontsize',8); 
end
text(f4.ax(5),0.65,0.3,{'in bound 2'},'unit','normalized');
xlim(f4.ax(5),[-1 6]);
hold(f4.ax(5),'off'); 

hold(f4.ax(6),'on');
text(f4.ax(6),0.65,0.4,strcat({'evt rate: '},sprintf('%.1f',rtouttmr2)),'unit',...
     'normalized');
for i = 1: length(h6.Values)
    text(f4.ax(6),h6.BinEdges(i),h6.Values(i)+f4.ax(6).YLim(2)/30,sprintf('%.1f',rtouttmr2m(i)),...
         'fontsize',8); 
end
xlim(f4.ax(6),[-1 6]);
hold(f4.ax(6),'off'); 
        
supertit(f4.ax,prefix);        
        
% save figure
print(f4.fig,'-dpdf',strcat(figpath,'/evt_freq_mag_diff_bounds.',prefix,'.dep',...
      num2str(depran),'.pdf'));        
        
       

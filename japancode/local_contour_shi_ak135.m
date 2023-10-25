% function local_contour_shi_ak135
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
depran = 15;
% depran = inf;

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
    evtall = Jmacatlogread(regflag,depran,slab);

else    % i.e. no need to recalculate, load the existing results
    % load the existing data
    data = load(strcat(evtpath,'/regevt',prefix,'dep',num2str(depran),'.mat'));
    evtall = data.evtall;
end
% ignore the ones with magnitude undetermined
evtall = evtall(evtall(:,11)~=1e6, :);

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
evtall = evtall(evtall(:,1) >= obara(1,1) & evtall(:,1) <= obara(end,1), :);

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
        132.6610   33.5583;
        132.8415   33.4803;
        132.8503   33.2645;
        132.7540   32.8785;
        132.2677   32.6653;
        132.1118   32.6378;
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
csp1x = linspace(134.1,133.7,1000);
csp1y = linspace(33.6,34.9,1000);
csp1 = [csp1x' csp1y'];

csp2x = linspace(132.6,131.7,1000);
csp2y = linspace(33.0,33.8,1000);
csp2 = [csp2x' csp2y'];

% assume tremors are on the slab interface, interpolate to get their depth
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
obara(:,8) = F(obara(:,6),obara(:,7));

%%%% for bound 1 events and their cross-section profile
% define the origin as the lower right point
lon0 = csp1(1,1);
lat0 = csp1(1,2);

% convert the cross-section profile to relative distance 
[csp1r(:,1),csp1r(:,2)] = absloc2relaloc(csp1(:,1),csp1(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(csp1r(end,2),csp1r(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
csp1r(:,3) = F(csp1(:,1),csp1(:,2));
% rotate the coordinates to profile axis 
[csp1r(:,4),csp1r(:,5)] = coordinate_rot(csp1r(:,1),csp1r(:,2),rotang,0,0);

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
lon0 = csp2(1,1);
lat0 = csp2(1,2);

% convert the cross-section profile to relative distance 
[csp2r(:,1),csp2r(:,2)] = absloc2relaloc(csp2(:,1),csp2(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(csp2r(end,2),csp2r(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
csp2r(:,3) = F(csp2(:,1),csp2(:,2));
% rotate the coordinates to profile axis 
[csp2r(:,4),csp2r(:,5)] = coordinate_rot(csp2r(:,1),csp2r(:,2),rotang,0,0);

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

% choose one
% velmod = jma2001;
% velmod = ukawa83;
velmod = 'ak135';

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


%% contour for tremor azimuth and ray p 
i=1;
raypsta1 = rayptmr1(:,i)*pi/180;
azsta1 = aztmr1(:,i);
lat = 32.3:0.01:35.2;
lon = 131:0.01:135;
[longrd, latgrd] = meshgrid(lon,lat);
raypgrd = griddata(tmrsel(:,6),tmrsel(:,7),raypsta1,longrd,latgrd, 'cubic');
azgrd = griddata(tmrsel(:,6),tmrsel(:,7),azsta1,longrd,latgrd, 'cubic');

f.fig=figure;
f.fig.Renderer='Painters';
widin = 15;  % maximum width allowed is 8.5 inches
htin = 12;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f.ax = gca;
hold(f.ax,'on');
grid(f.ax, 'on');
f.ax.Box = 'on';
title(f.ax,strcat(prefix,{'-tremor-'},velflag),'fontsize',12);
if regflag == 1  % means western shikoku
    plot(f.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
         [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
axis(f.ax, 'equal');
if regflag == 1
    xlim(f.ax,[131.8 134.6]);
    ylim(f.ax,[32.8 34.9]);
    xticks(f.ax,131.8:0.25:134.6);
    yticks(f.ax,32.8:0.25:34.9);
    f.ax.XMinorGrid = 'on';
    f.ax.XMinorTick = 'on';
    f.ax.YMinorGrid = 'on';
    f.ax.YMinorTick = 'on';
elseif regflag == 2
    xlim(f.ax,[134.5 137.2]);
    ylim(f.ax,[33 35.5]);
end
scatter(f.ax,evtbnd1(:,8),evtbnd1(:,9),2,[0.5 0.8 1],'filled','o');
scatter(f.ax,obara(:,6),obara(:,7),2,[0.7 0.7 0.7],'filled','o');

scatter(f.ax,stlosel1(i), stlasel1(i),60,'y','filled','^','markeredgec','k');
text(f.ax,stlosel1(i), stlasel1(i)+0.03, stnmsel1(i,:), 'fontsize',6);

[~,ct1]=contourf(f.ax,longrd,latgrd,raypgrd,floor(min(raypsta1)):0.2:ceil(max(raypsta1)),'k-',...
        'showtext','on','LabelSpacing',500); %,'linec','k'
caxis(f.ax,[min(raypsta1),max(raypsta1)]);
[~,ct2]=contour(f.ax,longrd,latgrd,azgrd,floor(min(azsta1)):2:ceil(max(azsta1)),'w--',...
        'linew',1.5);%,'linec','b','showtext','on','LabelSpacing',300
legend(f.ax,[ct1,ct2],{'ray p, interval 0.2 s/deg','azi, interval 2 deg'},'location',...
       'northwest','fontsize',12);
% save figure
print(f.fig,'-depsc2',strcat(figpath,'/tmr_raypazi_contour.',prefix,'.dep',...
      num2str(depran),'_',velflag,'.eps'));   


%% contour for events azimuth and ray p 
i=1;
raypsta1 = raypevt1(:,i)*pi/180;
azsta1 = azevt1(:,i);
lat = 32.3:0.01:35.2;
lon = 131:0.01:135;
[longrd, latgrd] = meshgrid(lon,lat);
raypgrd = griddata(evtbnd1(:,8),evtbnd1(:,9),raypsta1,longrd,latgrd, 'cubic');
azgrd = griddata(evtbnd1(:,8),evtbnd1(:,9),azsta1,longrd,latgrd, 'cubic');

f.fig=figure;
f.fig.Renderer='Painters';
widin = 15;  % maximum width allowed is 8.5 inches
htin = 12;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f.ax = gca;
hold(f.ax,'on');
grid(f.ax, 'on');
f.ax.Box = 'on';
title(f.ax,strcat(prefix,{'-event-'},velflag),'fontsize',12);
if regflag == 1  % means western shikoku
    plot(f.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
         [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
axis(f.ax, 'equal');
if regflag == 1
    xlim(f.ax,[133.5 134.7]);
    ylim(f.ax,[33.4 34.2]);
    xticks(f.ax,133.5:0.25:134.7);
    yticks(f.ax,33.4:0.25:34.2);
    f.ax.XMinorGrid = 'on';
    f.ax.XMinorTick = 'on';
    f.ax.YMinorGrid = 'on';
    f.ax.YMinorTick = 'on';
elseif regflag == 2
    xlim(f.ax,[134.5 137.2]);
    ylim(f.ax,[33 35.5]);
end
scatter(f.ax,evtbnd1(:,8),evtbnd1(:,9),2,[0.5 0.8 1],'filled','o');
scatter(f.ax,obara(:,6),obara(:,7),2,[0.7 0.7 0.7],'filled','o');

scatter(f.ax,stlosel1(i), stlasel1(i),60,'y','filled','^','markeredgec','k');
text(f.ax,stlosel1(i), stlasel1(i)+0.03, stnmsel1(i,:), 'fontsize',6);

[~,ct1]=contourf(f.ax,longrd,latgrd,raypgrd,floor(min(raypsta1)):0.2:ceil(max(raypsta1)),'k-',...
        'showtext','on','LabelSpacing',500); %,'linec','k'
caxis(f.ax,[min(rayptmr1(:,i)*pi/180),max(rayptmr1(:,i)*pi/180)]);
[~,ct2]=contour(f.ax,longrd,latgrd,azgrd,floor(min(azsta1)):2:ceil(max(azsta1)),'w--',...
        'linew',1.5);%,'linec','b','showtext','on','LabelSpacing',300
legend(f.ax,[ct1,ct2],{'ray p, spacing 0.2 s/deg','azi, spacing 2 deg'},'location',...
       'northwest','fontsize',12);
% save figure
print(f.fig,'-depsc2',strcat(figpath,'/evt_raypazi_contour.',prefix,'.dep',...
      num2str(depran),'_',velflag,'.eps')); 






















% function evt_piercing_point_kii.m
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
% First created date:   2020/03/17
% Last modified date:   2020/03/17
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
% obara(:,8) = griddata(slab(:,1),slab(:,2),-slab(:,3),obara(:,6),obara(:,7),'cubic');

% define the origin as the mid lower right point
lon0 = cspmid(1,1);
lat0 = cspmid(1,2);

% convert the cross-section profile to relative distance 
[cspmidr(:,1),cspmidr(:,2)] = absloc2relaloc(cspmid(:,1),cspmid(:,2),lon0,lat0);

% compute the rotation angle counter-clockwise from y to the profile
rotang = rad2deg((atan2(cspmidr(end,2),cspmidr(end,1)))) - 90;

% also interpolate the slab along the profile, because along the profile, 3D slab is only a line
cspmidr(:,3) = F(cspmid(:,1),cspmid(:,2));
% cspmidr(:,3) = griddata(slab(:,1),slab(:,2),-slab(:,3),cspmid(:,1),cspmid(:,2),'cubic');

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
stasel = [];
for i = 1: length(indsel)
    stasel(i).nm = stnm(indsel(i),:);
    stasel(i).lo = stlo(indsel(i));    % station longitude
    stasel(i).la = stla(indsel(i));    % station latitude
    stasel(i).el = stel(indsel(i));    % station elevation
end

% stnmsel = stnm(indsel,:);    % station name, string
% stlosel = stlo(indsel);    % station longitude
% stlasel = stla(indsel);    % station latitude
% stelsel = stel(indsel);    % station elevation


%% calculate the piercing points of events on the slab interface
recalflag = 0;
dslab = 0.002;
evtsel = evtbnd1;
kdense = 10;
velmod = 'jma2001';
% velmod = 'ukawa83';
% velmod = 'ak135';
rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% ttsevt =rayinfo.tts;
% raypevt = rayinfo.rayp;
% distevt = rayinfo.dist;
% pierptevt = rayinfo.pierpt ;
% poorresevt = rayinfo.poorres;
% edistevt = rayinfo.edist ;

ttsevt =rayinfo.ttsevt;
raypevt = rayinfo.raypevt;
distevt = rayinfo.distevt;
pierptevt = rayinfo.pierptevt;
poorresevt = rayinfo.poorres;
edistevt = rayinfo.edistevt;


%% choose a station to plot the piercing points
%        'N.KRTH';
%        'N.TKEH';
%        'N.KAWH';
%        'N.OWSH';
%        'N.HYSH';
%        'N.MGWH'
sta = ['N.INMH';
       'N.HRKH';
       'N.NKMH';
       'N.HNZH';
       'N.TKWH'];
% sta = 'N.INMH';
for ista = 1:size(sta,1)
   
for i = 1:size(stasel,2)
    if isequal(stasel(i).nm, sta(ista,:))
        stanum = i;
    end
end

nsta = size(pierptevt,1);
nevt = size(pierptevt,2);

pierpt(:,:) = pierptevt(stanum,:,:);


% %% plot the 3d view of piercing points and tremors 
% f.fig=figure;
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 8;   % maximum height allowed is 11 inches
% set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);
% 
% f.ax = gca; 
% f.ax = plt_tmr_evt_pier_point3D(f.ax,regflag,slab,obara,pierpt,stasel(stanum));
% 
% title(f.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% 
% % save figure
% savefig(f.fig,strcat(figpath,'/tmr_pierpt_3dview.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.fig'));
    
   

% %% obtain the statistics of distance between tremors and piercing points
% range = 2.5;      % max dist, [km]
% plotflag = 1;   % output the plots
% [ind,eucdist,ind2,eucdist2,numtmr,f1,f2,f3,~,f5]=...
%     tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,pierpt,stasel(stanum),range);
% 
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
% 
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f3.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f3.fig,'-dpdf',strcat(figpath,'/num_tmr_pierpt_onslab_inran.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% % title(f4.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% % print(f4.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_mindist_onslab_inran.',prefix,'_',stasel(stanum).nm,...
% %         '_',velmod,'.pdf'));    
% 
% title(f5.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f5.fig,'-dpdf',strcat(figpath,'/num_tmr_pierpt_inran_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));


    
% %% compare the diffence between different velocity models
% %%%% usinng function 'plt_tmr_peript_dist_diff.m'
% range = 2.5;
% plotflag = 1;   % do not output the plots
% 
% velmod = 'jma2001';
% rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% pierptevt = rayinfo.pierptevt;
% pierptjma(:,:) = pierptevt(stanum,:,:);
% [~,distjma,~,~,~,f1,f2,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                              pierptjma,stasel(stanum),range);
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
%     
% velmod = 'ukawa83';
% rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% pierptevt = rayinfo.pierpt;
% pierptuka(:,:) = pierptevt(stanum,:,:);
% [~,distuka,~,~,~,f1,f2,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                              pierptuka,stasel(stanum),range);                                         
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
%     
% velmod = 'ak135';                                         
% rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% pierptevt = rayinfo.pierpt;
% pierptak(:,:) = pierptevt(stanum,:,:);
% [~,distak,~,~,~,f1,f2,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                              pierptak,stasel(stanum),range);
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%  
%     
% diffuka = distuka-distjma;
% diffak = distak-distjma;
% [f1,f2]=plt_tmr_pierpt_dist_diff(regflag,obara,evtsel,stasel(stanum),diffuka);
% [f3,f4]=plt_tmr_pierpt_dist_diff(regflag,obara,evtsel,stasel(stanum),diffak);
% 
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ukawa83-jma2001'));
% print(f1.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist.',prefix,'_',stasel(stanum).nm,...
%         '_','ukawa83-jma2001','.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ukawa83-jma2001'));
% print(f2.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_','ukawa83-jma2001','.pdf'));
%     
% title(f3.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ak135-jma2001'));
% print(f3.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist.',prefix,'_',stasel(stanum).nm,...
%         '_','ak135-jma2001','.pdf'));
%     
% title(f4.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ak135-jma2001'));
% print(f4.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_','ak135-jma2001','.pdf'));



%% compare the diffence between different velocity models
%%%% usinng function 'plt_tmr_peript_diffmodel.m'
range = 2.5;
plotflag = 1;   % do not output the plots

velmod = 'jma2001';
rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
pierptevt = rayinfo.pierptevt;
pierptjma(:,:) = pierptevt(stanum,:,:);
% [~,distjma,~,~,f1,f2,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                              pierptjma,stasel(stanum),range);
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
    
    
velmod = 'ukawa83';
rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
pierptevt = rayinfo.pierpt;
pierptuka(:,:) = pierptevt(stanum,:,:);
% [~,distuka,~,~,f1,f2,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                              pierptuka,stasel(stanum),range);                                         
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
    
    
velmod = 'ak135';                                         
rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
pierptevt = rayinfo.pierpt;
pierptak(:,:) = pierptevt(stanum,:,:);
% [~,distak,~,~,f1,f2,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                              pierptak,stasel(stanum),range);
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f1.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_onslab.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% print(f2.fig,'-dpdf',strcat(figpath,'/tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
    
plotflag = 1;

[f1,f2,diffuka] = plt_tmr_pierpt_diffmodel(regflag,plotflag,pierptjma,pierptuka,obara,evtsel,stasel(stanum));

[f3,f4,diffak] = plt_tmr_pierpt_diffmodel(regflag,plotflag,pierptjma,pierptak,obara,evtsel,stasel(stanum));

perc95uka = prctile(diffuka,95);
perc95ak = prctile(diffak,95);

referror(ista) = 1/2*(perc95uka + perc95ak);

title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ukawa83-jma2001'));
print(f1.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist.',prefix,'_',stasel(stanum).nm,...
        '_','ukawa83-jma2001','.pdf'));
    
title(f2.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ukawa83-jma2001'));
print(f2.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
        '_','ukawa83-jma2001','.pdf'));
    
title(f3.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ak135-jma2001'));
print(f3.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist.',prefix,'_',stasel(stanum).nm,...
        '_','ak135-jma2001','.pdf'));
    
title(f4.ax,strcat(prefix,'---',stasel(stanum).nm,'---','ak135-jma2001'));
print(f4.fig,'-dpdf',strcat(figpath,'/diff_tmr_pierpt_dist_hist.',prefix,'_',stasel(stanum).nm,...
        '_','ak135-jma2001','.pdf'));

end
















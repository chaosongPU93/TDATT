% function evt_samepathas_tmr_kii.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to obtain the piercing point of the ray path from the 
% event to station on the slab interface. First, creat a kd-tree model 
% (KDtreeSearcher) based on the slab grid data (after interpolation). Second,
% use rangesearch to find all tremors that are 'close' enough to the piercing
% point of each event at a fixed station. This range partially comes from the
% difference between piercing points obtained with different velocity models.
% Because we want to minimize the effect from velocity models, thus stations
% close to tremors and events (in perfect case on top of the tremors) are
% chosen only. Let's say 95 percentile of the piercing points would have a
% location difference of 'dmod' km, then the searching range R should be no
% smaller than dmod, i.e., R>=dmod, here assuming 95% is significant enough.
% 
% 
%
% kdtree:
% 1.Kd-trees divide your data into nodes with at most BucketSize (default
%   is 50) points per node, based on coordinates (as opposed to categories).
%
%
%
% 
% range search:
% 1. Find all neighbors within specified distance using searcher object.
% 2. Basically the same syntax as knnsearch, Idx = rangesearch(Mdl,Y,r) searches
%   for all neighbors (i.e., points, rows, or observations) in Mdl.X within radius
%   r of each point (i.e., row or observation) in the query data Y using an
%   exhaustive search or a Kd-tree. rangesearch returns Idx, which is a column
%   vector of the indices of Mdl.X within r units.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/05/07
% Last modified date:   2020/05/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
clc;
clear;
close all;

% set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

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


%% read the piercing points of events on the slab interface
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


%% choose a station to obtain the distance between tremors and piercing points
sta = 'N.INMH';
% sta = 'N.HRKH';
% sta = 'N.HNZH';
% sta = 'N.TKWH';
% sta = 'N.KRTH';
% sta = 'N.KAWH';

for i = 1:size(stasel,2)
    if isequal(stasel(i).nm, sta)
        stanum = i;
    end
end

nsta = size(pierptevt,1);
nevt = size(pierptevt,2);

pierpt(:,:) = pierptevt(stanum,:,:);


% %% a generic analysis
% range = 0.5:0.5:10;
% fracin = zeros(size(range));
% N = zeros(size(range));
% h = [];
% for j = 1: length(range)
% 
%     plotflag = 0;   % do not output the plots
%     % do a range search instead of knn, ouput the ind and eucdist of tremors within the range round
%     % each piercing point 
%     [~,~,ind,eucdist,~,~,~,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
%                                 pierpt,stasel(stanum),range(j));
%                             
%     numtmr = zeros(size(eucdist));
%     tmrinind = cell(size(eucdist)); % cell stores the index of tremor have co-time events, can be empty
%     tmroutind = cell(size(eucdist));% cell stores the index of tremor have co-path events, can be empty
%     indcopath = [];     % stores the index of events that have co-path tremors
%     indcotime = [];     % stores the index of events that have co-path tremors occuring at the same time
%     for i = 1: size(eucdist,1)
%         % i=4;
%         numtmr(i) = size(eucdist{i},2);
%         if numtmr(i) ~= 0
%             indcopath = [indcopath; i];
%             tmrind = ind{i};
%             tmrind = tmrind';
%             % obtain the index of tremor that shares the same time as
%             % the target event (i.e., occurred 1 hr or less before event); and the tremor that
%             % are outside this time period
%             [indin, indout] = check_simult_evt_tmr(obara(1,1),obara,evtsel(i,:),tmrind);
%             if ~isempty(indin)
%                 indcotime = [indcotime; i];
%                 tmrinind{i} = indin';
%             else
%                 tmrinind{i} = [];
%             end
%             tmroutind{i} = indout';
%         else
%             tmrinind{i} = [];
%             tmroutind{i} = [];
%         end
%     end
% 
%     fracin(j) = size(indcotime,1)/size(indcopath,1)*100;
%     N(j) = size(indcopath,1);
%     evtintmr = evtsel(indcotime,:);       % event occur at the same time as tremor
%     h1 = sum(evtintmr(:,11)>=1);
%     h2 = sum(evtintmr(:,11)<1);
%     h = [h;h1 h2];
% 
% end


% %% get the difference between velocity models
% velmod = 'jma2001';
% rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% pierptevt = rayinfo.pierptevt;
% pierptjma(:,:) = pierptevt(stanum,:,:);
%        
% velmod = 'ukawa83';
% rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% pierptevt = rayinfo.pierpt;
% pierptuka(:,:) = pierptevt(stanum,:,:);
% 
% velmod = 'ak135';                                         
% rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod);
% pierptevt = rayinfo.pierpt;
% pierptak(:,:) = pierptevt(stanum,:,:);
% 
% plotflag = 1;
% [~,~,diffuka] = plt_tmr_pierpt_diffmodel(regflag,plotflag,pierptjma,pierptuka,obara,evtsel,...
%                                            stasel(stanum));
% 
% [~,~,diffak] = plt_tmr_pierpt_diffmodel(regflag,plotflag,pierptjma,pierptak,obara,evtsel,...
%                                           stasel(stanum));
% perc95uka = prctile(diffuka,95);
% perc95ak = prctile(diffak,95);
% 
% referror = max(perc95uka + perc95ak);


% %% plot the variation
% velmod = 'jma2001';
% 
% f1.fig=figure;
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 7;   % maximum height allowed is 11 inches
% set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% f1.ax = gca;
% hold(f1.ax,'on');
% f1.ax.Box = 'on';
% 
% yyaxis(f1.ax,'right');
% H = bar(f1.ax,range,h,1,'stacked','facea',0.7);
% H(1).FaceColor = [0.3 0.3 0.3];
% H(2).FaceColor = [0.7 0.7 0.7];
% H(1).EdgeColor = 'w';
% H(2).EdgeColor = 'w';
% ylabel(f1.ax,'Num of evts same time as tmr w/ search range');
% 
% yyaxis(f1.ax,'left');
% plot(f1.ax,range,fracin,'.-','linew',1,'color','b','markers',16);
% for i = 1: length(range)
%     text(f1.ax,range(i)-0.25,fracin(i)+1e-2,num2str(N(i)),'fontsize',9);
% end
% xlabel(f1.ax,'Radius of search range (km)');
% ylabel(f1.ax,'Frac of evts same time as tmr w/ search range (%)');
% 
% f1.ax.YAxis(1).Color = [0.15 0.15 0.15];
% f1.ax.YAxis(2).Color = [0.15 0.15 0.15];
% f1.ax.XAxis.Limits = [0 10.5];
% ylim1 = f1.ax.YAxis(1).Limits;
% f1.ax.YAxis(1).Limits(2) = 0.2*(ylim1(2)-ylim1(1))+ylim1(2);
% ylim2 = f1.ax.YAxis(2).Limits;
% f1.ax.YAxis(2).Limits(2) = 0.2*(ylim2(2)-ylim2(1))+ylim2(2);
% 
% plot(f1.ax,[referror, referror], [0,20], 'k--','linew',1.5);
% text(f1.ax,referror+0.25, 0.8*f1.ax.YAxis(1).Limits(2), 'max. 95% vel. difference');
% 
% legend(f1.ax,H,{'Mag>=1','Mag<1'},'location','northwest');
% title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% 
% print(f1.fig,'-dpdf',strcat(figpath,'/evt_samepath_tmr_distvari.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
    

%% choose a proper range from above and another constraint
% this range also comes from the statistical result of 'diffdist' otuput from
% 'plt_tmr_pierpt_diffmodel.m'
distran = 4;
plotflag = 0;   % do not output the plots

% do a range search instead of knn, ouput the ind and eucdist of tremors within the range round
% each piercing point
[~,~,ind,eucdist,~,kdtreetmr,~,~,~,~]=tmr_pierpt_dist_onslab(regflag,plotflag,slab,obara,evtsel,...
                                pierpt,stasel(stanum),distran);
                            
numtmr = zeros(size(eucdist));
tmrinind = cell(size(eucdist));% cell stores the index of tremor have co-time events, can be empty
tmroutind = cell(size(eucdist));% cell stores the index of tremor have co-path events, can be empty
indcopath = [];     % stores the index of events that have co-path tremors
indcotime = [];     % stores the index of events that have co-path tremors occuring at the same time 
for i = 1: size(eucdist,1)
% i=4;
    numtmr(i) = size(eucdist{i},2);
    if numtmr(i) ~= 0
        indcopath = [indcopath; i];
        tmrind = ind{i};
        tmrind = tmrind';
        % obtain the index of tremor that shares the same time as
        % the target event (i.e., occurred 1 hr or less before event); and the tremor that
        % are outside this time period
        [indin, indout] = check_simult_evt_tmr(obara(1,1),obara,evtsel(i,:),tmrind);
        if ~isempty(indin)
            indcotime = [indcotime; i];
            tmrinind{i} = indin';
        else
            tmrinind{i} = [];
        end
        tmroutind{i} = indout';
    else
        tmrinind{i} = [];
        tmroutind{i} = []; 
    end
end

evtcopath = evtsel(indcopath,:);    % event have co-path tremors
evtintmr = evtsel(indcotime,:);       % event occur at the same time as tremor
indncotime = setdiff(indcopath,indcotime);
evtouttmr = evtsel(indncotime,:);      % event not occur at the same time as tremor
tmrinevt = cell(size(indcotime));   % cell array 
for i = 1:length(indcotime)
    indtemp = tmrinind{indcotime(i)};
    tmrinevt{i} = obara(indtemp,:);
end

fracin = size(evtintmr,1)/size(evtcopath,1)*100;


%% Find the tremors within +-24 hours of the target events
% to see the trend of tremor migration if there is a trend of tremor migration, also to make sure
% that the associated co-time tremor is not an entirely isolated

% relative time in days for the first day of tremor catalog           
d1 = datetime(obara(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

% relative origin time in days of all tremors
ditmr = datetime(obara(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
ottmr = caldays(between(d1,ditmr,'days')) + obara(:,5)/24;      % unit is day

tmrbkgd = cell(size(indcotime));   % cell array
for i = 1:length(indcotime)
    diobj = datetime(evtintmr(i,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    otobj = caldays(between(d1,diobj,'days')) + (evtintmr(i,5)+evtintmr(i,6)/60+...
                    evtintmr(i,7)/3600)/24;
    tmrbkgd{i} = obara(ottmr>=otobj-1 & ottmr<=otobj+1, :);        
end


%% Find the other events that within the same radius, ie same path, but different time as target 
% this is better to be done with a new kdtree created by all piercing points
if regflag == 1
    prefix = 'shikoku';
    lat = 32.5:0.01:35.5;
    lon = 131.5:0.01:135;
elseif regflag == 2
    prefix = 'kii';
    lat = 32.5:0.01:35;
    lon = 134:0.01:137.5;
end
lat0 = 0.5*(lat(1)+lat(end));
lon0 = 0.5*(lon(1)+lon(end));
% relative coordinates of piercing points
[pierptr(:,1),pierptr(:,2)] = absloc2relaloc(pierpt(:,1),pierpt(:,2),lon0,lat0);
pierptr(:,3) = pierpt(:,3);

% create a kd tree for all piercing points
kdtreeppt = KDTreeSearcher(pierptr,'Distance','euclidean','BucketSize',50);

% rangesearch, find the all other piercing points within the dist range to each target piercing point
% i.e. within this radius, we think their ray paths are the same
% the result is a cell array
ppttarget = pierptr(indcotime, :);
[ind3,eucdist3] = rangesearch(kdtreeppt,ppttarget,distran);

%%% relative time in days for the first day of tremor catalog           
d1 = datetime(obara(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

numevt = zeros(size(eucdist));  % number of copath events each target event has
evtbef = cell(size(eucdist3));  % copath events occurred 12 hr before the target event    
evteq = cell(size(eucdist3));   % copath events occurred +-12 hr of the target event 
evtaft = cell(size(eucdist3));  % copath events occurred 12 hr after the target event
evtbefind = cell(size(eucdist3));  % copath events occurred 12 hr before the target event    
evteqind = cell(size(eucdist3));   % copath events occurred +-12 hr of the target event 
evtaftind = cell(size(eucdist3));  % copath events occurred 12 hr after the target event
indcophevt = [];     % stores the index of target events that have other copath events

for i = 1: size(eucdist3,1)
% i=4;
    numevt(i) = size(eucdist3{i},2);
    if numevt(i) ~= 0
        indcophevt = [indcophevt; i];
        evtind = ind3{i};
        evtind = evtind';
        evtcophevt = evtsel(evtind,:);
        
        % relative origin time in days of all copath evts
        dievt = datetime(evtcophevt(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
        otevt = caldays(between(d1,dievt,'days')) + (evtcophevt(:,5)+evtcophevt(:,6)/60+...
                evtcophevt(:,7)/3600)/24;      % unit is day
        
        % relative origin time in days of all target events
        diobj = datetime(evtintmr(i,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
        otobj = caldays(between(d1,diobj,'days')) + (evtintmr(i,5)+evtintmr(i,6)/60+...
                evtintmr(i,7)/3600)/24;
        
        evtbefind{i} = evtind(otevt<otobj-12/24);
        evtbef{i} = evtcophevt(otevt<otobj-12/24, :);
        if size(evtbef{i},1) == 0
            evtbef{i} = [];
            evtbefind{i} = [];
        end
        
        evteqind{i} = evtind(otevt>=otobj-12/24 & otevt<=otobj+12/24 & otevt~=otobj);
        evteq{i} = evtcophevt(otevt>=otobj-12/24 & otevt<=otobj+12/24 & otevt~=otobj, :);
        if size(evteq{i},1) == 0
            evteq{i} = [];
            evteqind{i} = [];
        end
        
        evtaftind{i} = evtind(otevt>otobj+12/24);
        evtaft{i} = evtcophevt(otevt>otobj+12/24, :);
        if size(evtaft{i},1) == 0
            evtaft{i} = [];
            evtaftind{i} = [];
        end
    else
        evtbef{i} = [];
        evteq{i} = [];
        evtaft{i} = [];
        evtbefind{i} = [];
        evteqind{i} = [];
        evtaftind{i} = [];
    end
end

%%
% for these evtbef or evtaft, i also want to distinguish those events that the closest tremor before
% it within a distance range; if some tremor occurs outside this distance range, we thnik whenever
% it occurs, it would affect this event; but if within the distance, we want to distinguish the time
% of closest tremor is more than 5 hr or less 5 hr? Pressumably it is the safest if the cloest
% tremor occurred more than 5 hr or this event.

distran = 4;  % km, here i choose the same range to keep consistency
timeran = 5;  % hr, the time range to divide the closest tremors before query event in time

evtbefgt = cell(size(evtbef));
evtbefle = cell(size(evtbef));
for i = 1: size(evtbef,1)
    tmpevt = evtbef{i};
    [~,tmpind,~] = intersect(evtsel,tmpevt,'row');
    
    tmp1 = [];
    tmp2 = [];
    for j = 1: size(tmpevt,1)
        % find the all the tremors that occurred within the range of the piercing point of each evtbef 
        [ind4,~] = rangesearch(kdtreetmr,pierptr(tmpind(j),:),distran);  
        tmptmr = obara(ind4{1},:);
        
        % difference in time of the closest tremor before the query event
        [objtmr, mindt] = class_tmr_bef_evt(obara(1,1),tmptmr,tmpevt(j,:));
        if ~isempty(mindt)  
            if mindt > timeran
                tmp1 = [tmp1; tmpevt(j,:)];
            else
                tmp2 = [tmp2; tmpevt(j,:)];
            end
        else    % meaning all tremors occurs later than this event, so it is also safe
            tmp1 = [tmp1; tmpevt(j,:)];
        end
    end
    evtbefgt{i} = tmp1;
    evtbefle{i} = tmp2;
end

evtaftgt = cell(size(evtaft));
evtaftle = cell(size(evtaft));
for i = 1: size(evtaft,1)
    tmpevt = evtaft{i};
    [~,tmpind,~] = intersect(evtsel,tmpevt,'row');
    
    tmp1 = [];
    tmp2 = [];
    for j = 1: size(tmpevt,1)
        % find the all the tremors that occurred within the range of the piercing point of each evtbef 
        [ind4,~] = rangesearch(kdtreetmr,pierptr(tmpind(j),:),distran);  
        tmptmr = obara(ind4{1},:);
        
        % difference in time of the closest tremor before the query event
        [objtmr, mindt] = class_tmr_bef_evt(obara(1,1),tmptmr,tmpevt(j,:));
        if ~isempty(mindt)  
            if mindt > timeran
                tmp1 = [tmp1; tmpevt(j,:)];
            else
                tmp2 = [tmp2; tmpevt(j,:)];
            end
        else    % meaning all tremors occurs later than this event, so it is also safe
            tmp1 = [tmp1; tmpevt(j,:)];
        end
    end
    evtaftgt{i} = tmp1;
    evtaftle{i} = tmp2;
end

% try another distance range
distran = 10;  % km, here i choose the same range to keep consistency
timeran = 5;  % hr, the time range to divide the closest tremors before query event in time

evtbefgt10 = cell(size(evtbef));
evtbefle10 = cell(size(evtbef));
for i = 1: size(evtbef,1)
    tmpevt = evtbef{i};
    [~,tmpind,~] = intersect(evtsel,tmpevt,'row');
    
    tmp1 = [];
    tmp2 = [];
    for j = 1: size(tmpevt,1)
        % find the all the tremors that occurred within the range of the piercing point of each evtbef 
        [ind4,~] = rangesearch(kdtreetmr,pierptr(tmpind(j),:),distran);  
        tmptmr = obara(ind4{1},:);
        
        % difference in time of the closest tremor before the query event
        [objtmr, mindt] = class_tmr_bef_evt(obara(1,1),tmptmr,tmpevt(j,:));
        
        if ~isempty(mindt)  
            if mindt > timeran
                tmp1 = [tmp1; tmpevt(j,:)];
            else
                tmp2 = [tmp2; tmpevt(j,:)];
            end
        else    % meaning all tremors occurs later than this event, so it is also safe
            tmp1 = [tmp1; tmpevt(j,:)];
        end
    end
    evtbefgt10{i} = tmp1;
    evtbefle10{i} = tmp2;
end

evtaftgt10 = cell(size(evtaft));
evtaftle10 = cell(size(evtaft));
for i = 1: size(evtaft,1)
    tmpevt = evtaft{i};
    [~,tmpind,~] = intersect(evtsel,tmpevt,'row');
    
    tmp1 = [];
    tmp2 = [];
    for j = 1: size(tmpevt,1)
        % find the all the tremors that occurred within the range of the piercing point of each evtbef 
        [ind4,~] = rangesearch(kdtreetmr,pierptr(tmpind(j),:),distran);  
        tmptmr = obara(ind4{1},:);
        
        % difference in time of the closest tremor before the query event
        [objtmr, mindt] = class_tmr_bef_evt(obara(1,1),tmptmr,tmpevt(j,:));
        if ~isempty(mindt)  
            if mindt > timeran
                tmp1 = [tmp1; tmpevt(j,:)];
            else
                tmp2 = [tmp2; tmpevt(j,:)];
            end
        else    % meaning all tremors occurs later than this event, so it is also safe
            tmp1 = [tmp1; tmpevt(j,:)];
        end
    end
    evtaftgt10{i} = tmp1;
    evtaftle10{i} = tmp2;
end


%% A map for all target co-path and co-time tremor and event (piercing point) pairs
f1.fig=figure;
f1.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f1.ax = gca;
hold(f1.ax,'on');
f1.ax.Box = 'on';
grid(f1.ax, 'on');
if regflag == 1  % means western shikoku
    plot(f1.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
        [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f1.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
plot(f1.ax,[135.62,136.22],[34.78,33.88],'color',[0.6 0.6 0.6],'linew',2);
axis(f1.ax, 'equal');
if regflag == 1
    xlim(f1.ax,[131.5 133.5]);
    ylim(f1.ax,[32.5 34]);
    xticks(f1.ax,131.8:0.25:134.6);
    yticks(f1.ax,32.8:0.25:34.9);
    f1.ax.XMinorGrid = 'on';
    f1.ax.XMinorTick = 'on';
    f1.ax.YMinorGrid = 'on';
    f1.ax.YMinorTick = 'on';
elseif regflag == 2
    xlim(f1.ax,[134.5 136.4]);
    ylim(f1.ax,[33.1 34.8]);
    xticks(f1.ax,134.5:0.25:137);
    yticks(f1.ax,33.1:0.25:35.4);
    f1.ax.XMinorGrid = 'on';
    f1.ax.XMinorTick = 'on';
    f1.ax.YMinorGrid = 'on';
    f1.ax.YMinorTick = 'on';
end

lat = 31:0.01:36;
lon = 131:0.01:138;
[longrd, latgrd] = meshgrid(lon,lat);
depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');
contour(f1.ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.7 0.7 0.7]);

scatter(f1.ax,pierpt(:,1),pierpt(:,2),2,[0 1 1],'filled','s');% colored by dist from closest tremor

scatter(f1.ax,obara(:,6),obara(:,7),2,[0.8 0.8 0.8],'filled','o');

% scatter(f1.ax,evtintmr(:,8),evtintmr(:,9),15,'b','filled','s');

for i = 1: size(tmrinevt,1)
    tmrobj = tmrinevt{i};
    scatter(f1.ax,tmrobj(:,6),tmrobj(:,7),15,[0.4,0.4,0.4],'filled','o');
    scatter(f1.ax,pierpt(indcotime(i), 1),pierpt(indcotime(i), 2),15,'b','filled','s');
    plot(f1.ax,[tmrobj(:,6) pierpt(indcotime(i), 1)], [tmrobj(:,7) pierpt(indcotime(i), 2)], ...
         'k-');
end
scatter(f1.ax,stasel(stanum).lo, stasel(stanum).la,80,'y','filled','^','markeredgec','k');

text(f1.ax,stasel(stanum).lo+0.1, stasel(stanum).la-0.1, stasel(stanum).nm, 'fontsize',12, ...
     'color','k');
title(f1.ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));

print(f1.fig,'-dpdf',strcat(figpath,'/all_co-path-time_tmr_pierpt.',prefix,'_',stasel(stanum).nm,...
        '_',velmod,'.pdf'));
    

%% Separate maps for each target event (piercing point)
% subplot 1: co-path-time event and tremor, and tremor activity +-24 hours of the target events
% subplot 2: co-path-time event and tremor, and other co-path events but outside the tremor period
for i = 1: size(evtintmr,1)
    f2.fig=figure;
    f2.fig.Renderer='Painters';
    widin = 15;  % maximum width allowed is 8.5 inches
    htin = 7;   % maximum height allowed is 11 inches
    set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    
    f2.ax(1) = subplot(121);
    f2.ax(2) = subplot(122);
    
    % reposition
    set(f2.ax(1), 'position', [ 0.08, 0.08, 0.4, 0.9]);
    set(f2.ax(2), 'position', [ 0.55, 0.1, 0.4, 0.9]);
    
    % subplot 1: co-path-time event and tremor, and tremor activity +-24 hours of the target events
    ax = f2.ax(1);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    if regflag == 1  % means western shikoku
        plot(ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    elseif regflag == 2  % means kii pen
        plot(ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    end
    plot(ax,[135.62,136.22],[34.78,33.88],'color',[0.6 0.6 0.6],'linew',2);
    axis(ax, 'equal');
    if regflag == 1
        xlim(ax,[131.5 133.5]);
        ylim(ax,[32.5 34]);
    elseif regflag == 2
        xlim(ax,[134.5 136.4]);
        ylim(ax,[33.1 34.8]);
    end
    
    contour(ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.7 0.7 0.7]);
    
    scatter(ax,pierpt(:,1),pierpt(:,2),2,[0 1 1],'filled','s');% colored by dist from closest tremor
    
    scatter(ax,obara(:,6),obara(:,7),2,[0.8 0.8 0.8],'filled','o');

    diobj = datetime(evtintmr(i,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    otobj = caldays(between(d1,diobj,'days')) + (evtintmr(i,5)+evtintmr(i,6)/60+...
                    evtintmr(i,7)/3600)/24;
    
    tmrobj = tmrbkgd{i};
    ditmr = datetime(tmrobj(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
    ottmr = caldays(between(d1,ditmr,'days')) + tmrobj(:,5)/24;      % unit is day
        
    dt =  (ottmr-otobj)*24;
    scatter(ax,tmrobj(:,6),tmrobj(:,7),10,dt,'filled','o');
    colormap(ax,'jet');
    c=colorbar(ax);
    c.Label.String = strcat('Time of tremor wrt. event');
    c.Label.FontSize = 12;
    if ~isempty(ottmr)
        caxis(ax,[min(dt)-1 max(dt)+1]);
    end
        
    
    tmrobj = tmrinevt{i};
    scatter(ax,tmrobj(:,6),tmrobj(:,7),20,'k','filled','o');
    scatter(ax,pierpt(indcotime(i), 1),pierpt(indcotime(i), 2),25,'k','filled','s');
    plot(ax,[tmrobj(:,6) pierpt(indcotime(i), 1)], [tmrobj(:,7) pierpt(indcotime(i), 2)], ...
         'k-');
    
    scatter(ax,stasel(stanum).lo, stasel(stanum).la,80,'y','filled','^','markeredgec','k');
    
    text(ax,stasel(stanum).lo+0.1, stasel(stanum).la-0.1, stasel(stanum).nm, 'fontsize',12, ...
        'color','k');
    evtlabel = strcat(num2str(evtintmr(i,1)),num2str(evtintmr(i,5)), num2str(evtintmr(i,6)), ...
                      num2str(evtintmr(i,7)) );
                  
                  
    % subplot 2: co-path-time event and tremor, and other co-path events but outside the tremor period
    ax = f2.ax(2);
    hold(ax,'on');
    ax.Box = 'on';
    grid(ax, 'on');
    if regflag == 1  % means western shikoku
        plot(ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    elseif regflag == 2  % means kii pen
        plot(ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    end
    plot(ax,[135.62,136.22],[34.78,33.88],'color',[0.6 0.6 0.6],'linew',2);
    axis(ax, 'equal');
    if regflag == 1
        xlim(ax,[131.5 133.5]);
        ylim(ax,[32.5 34]);
    elseif regflag == 2
        xlim(ax,[134.5 136.4]);
        ylim(ax,[33.1 34.8]);
    end
    
    contour(ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.7 0.7 0.7]);
    
    scatter(ax,pierpt(:,1),pierpt(:,2),2,[0 1 1],'filled','s');% colored by dist from closest tremor
    
    scatter(ax,obara(:,6),obara(:,7),2,[0.8 0.8 0.8],'filled','o');
    
    tmppierpt = pierpt(evtbefind{i}, :);        
    scatter(ax,tmppierpt(:,1),tmppierpt(:,2),2,'b','filled','s');
    
    tmppierpt = pierpt(evtaftind{i}, :);
    scatter(ax,tmppierpt(:,1),tmppierpt(:,2),2,'r','filled','s');
    
    tmrobj = tmrinevt{i};
    scatter(ax,tmrobj(:,6),tmrobj(:,7),20,'k','filled','o');
    scatter(ax,pierpt(indcotime(i), 1),pierpt(indcotime(i), 2),25,'k','filled','s');
    plot(ax,[tmrobj(:,6) pierpt(indcotime(i), 1)], [tmrobj(:,7) pierpt(indcotime(i), 2)], ...
         'k-');
     
    supertit(f2.fig,strcat(prefix,'---',stasel(stanum).nm,'---',velmod,'---',evtlabel));
    
    print(f2.fig,'-depsc2',strcat(figpath,'/co-path-time_tmr_pierpt_bkgdtmr.',prefix,'_',stasel(stanum).nm,...
        '_',velmod,'_',evtlabel,'.eps'));

end


%% analyse the above plots to determine the indexes of copath-time events to be saved
if strcmp(sta, 'N.INMH')
    indsave = [2,3,4,7,8,9,10,11];  % 8
elseif strcmp(sta, 'N.HRKH')
    indsave = [3,5,7,8,9];  % 5
elseif strcmp(sta, 'N.HNZH')
    indsave = [5,6,9,10,12,14]; % 6
elseif strcmp(sta, 'N.TKWH')
    indsave = [5,7,9,10,11,13,14,18,19,20]; % 10
elseif strcmp(sta, 'N.KRTH')
    indsave = [2,3,4,5,6,7];    % 6
elseif strcmp(sta, 'N.KAWH')
    indsave = [4,5,7,9];    % 4
end    

% save the new matrixes
evtintmr = evtintmr(indsave, :);
tmrinevt = tmrinevt(indsave, :);
tmrbkgd = tmrbkgd(indsave, :);
evtbef = evtbef(indsave, :);
evtaft = evtaft(indsave, :);
evtbefgt = evtbefgt(indsave, :);
evtaftgt = evtaftgt(indsave, :);
evtbefle = evtbefle(indsave, :);
evtaftle = evtaftle(indsave, :);
evtbefgt10 = evtbefgt10(indsave, :);
evtaftgt10 = evtaftgt10(indsave, :);
evtbefle10 = evtbefle10(indsave, :);
evtaftle10 = evtaftle10(indsave, :);


% save to mat file
save(strcat(datapath,'/co-path-time_tmr_pierpt_',prefix,'_',stasel(stanum).nm,...
     '_',velmod,'.mat'),'evtintmr','tmrinevt','tmrbkgd','evtbef','evtaft','evtbefgt','evtaftgt',...
     'evtbefgt10','evtaftgt10','evtbefle','evtaftle','evtbefle10','evtaftle10');


% %% output the found these events
% % only the time precise to min is needed, ie. column 1 and 5,6
% 
% evtintmrtime = zeros(size(evtintmr,1),1);
% for i = 1: size(evtintmr,1)
%     if evtintmr(i,5) <= 9
%         hrstr = strcat('0',num2str(evtintmr(i,5)));
%     else
%         hrstr = num2str(evtintmr(i,5));
%     end
%     
%     if evtintmr(i,6) <= 9
%         minstr = strcat('0',num2str(evtintmr(i,6)));
%     else
%         minstr = num2str(evtintmr(i,6));
%     end
%     
%     timestr = strcat(num2str(evtintmr(i,1)),hrstr,minstr);
%     evtintmrtime(i) = str2num(timestr);
%     
% end
% 
% fid = fopen(strcat(evtpath,'/evts_samepath_timeas_tmr.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod),'w+');
% fprintf(fid,'%d %4d %2d %2d %2d %2d %5.2f %.4f %.4f %.2f %4.1f %7.1f %.4f %.4f \n',evtintmr');
% fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Below is not suggested, it is WRONG
%%%
% %% a generic analysis
% distsamp = 0.5:0.5:5;
% fracin = zeros(size(distsamp));
% N = zeros(size(distsamp));
% h = [];
% for i = 1: length(distsamp)
%     % keep in mind that here each tremor and event are paired
%     evtobj = evtsel(eucdist<=distsamp(i),:);    
%     tmrobj = tmrall(eucdist<=distsamp(i),:);
%     [evtin, ~] = class_simult_evt_tmr(obara(1,1),tmrobj,evtobj);
%     fracin(i) = size(evtin,1)/size(evtobj,1)*100;
%     N(i) = size(evtobj,1);
%     h1 = sum(evtin(:,11)>=1);
%     h2 = sum(evtin(:,11)<1);
%     h = [h;h1 h2];
% 
% end
% 
% %% plot the variation
% f2.fig=figure;
% ax = gca;
% hold(ax,'on');
% ax.Box = 'on';
% 
% yyaxis(ax,'right');
% H = bar(ax,distsamp,h,1,'stacked','facea',0.7);
% H(1).FaceColor = [0.3 0.3 0.3];
% H(2).FaceColor = [0.7 0.7 0.7];
% H(1).EdgeColor = 'w';
% H(2).EdgeColor = 'w';
% ylabel(ax,'Num of evts during tmr time within threshold');
% 
% yyaxis(ax,'left');
% plot(ax,distsamp,fracin,'.-','linew',1,'color','b','markers',16);
% for i = 1: length(distsamp)
%     text(ax,distsamp(i)-0.25,fracin(i)+1e-2,num2str(N(i)),'fontsize',9);
% end
% xlabel(ax,'Threshold of min dist (km)');
% ylabel(ax,'Frac of evts during tmr time within threshold (%)');
% 
% legend(ax,H,{'Mag>=1','Mag<1'},'location','northeast');
% 
% ax.YAxis(1).Color = [0.15 0.15 0.15];
% ax.YAxis(2).Color = [0.15 0.15 0.15];
% ax.XAxis.Limits = [0 5.5];
% ylim1 = ax.YAxis(1).Limits;
% ax.YAxis(1).Limits(2) = 0.2*(ylim1(2)-ylim1(1))+ylim1(2);
% ylim2 = ax.YAxis(2).Limits;
% ax.YAxis(2).Limits(2) = 0.2*(ylim2(2)-ylim2(1))+ylim2(2);
% title(ax,strcat(prefix,'---',stasel(stanum).nm,'---',velmod));
% 
% print(f2.fig,'-dpdf',strcat(figpath,'/evt_samepath_tmr_distvari.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod,'.pdf'));
%     
% 
% %% choose a distance threshold from above
% % keep in mind that here each tremor and event are paired
% 
% % events whose piercing points are as close as 1 km to tremors
% evtobj1 = evtsel(eucdist<=2,:);
% tmrobj1 = tmrall(eucdist<=2,:);
% 
% [evtin1, evtout1] = class_simult_evt_tmr(obara(1,1),tmrobj1,evtobj1);
% fracin1 = size(evtin1,1)/size(evtobj1,1);
% fracout1 = 1-fracin1;
% 
% 
% %% output the found these events
% % % only the time precise to min is needed, ie. column 1 and 5,6
% % 
% % evtin1time = zeros(size(evtin1,1),1);
% % for i = 1: size(evtin1,1)
% %     if evtin1(i,5) <= 9
% %         hrstr = strcat('0',num2str(evtin1(i,5)));
% %     else
% %         hrstr = num2str(evtin1(i,5));
% %     end
% %     
% %     if evtin1(i,6) <= 9
% %         minstr = strcat('0',num2str(evtin1(i,6)));
% %     else
% %         minstr = num2str(evtin1(i,6));
% %     end
% %     
% %     timestr = strcat(num2str(evtin1(i,1)),hrstr,minstr);
% %     evtin1time(i) = str2num(timestr);
% %     
% % end
% 
% fid = fopen(strcat(evtpath,'/evts_samepath_timeas_tmr.',prefix,'_',stasel(stanum).nm,...
%         '_',velmod),'w+');
% fprintf(fid,'%d %4d %2d %2d %2d %2d %5.2f %.4f %.4f %.2f %4.1f %7.1f %.4f %.4f \n',evtin1');
% fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















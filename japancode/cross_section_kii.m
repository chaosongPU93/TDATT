% function cross_section_kii.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to plot the cross-section of the earthquake distribution for Kii region
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
% First created date:   2020/03/17
% Last modified date:   2020/03/17
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

figure
% group scatter, scatter the same data set by group properties, without knowing how many options in
% that group variable. And this group variable can be any data type, char, string, number
gscatter(evtall(:,8),evtall(:,9),hypdet);


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

fid = fopen(strcat('/home/data2/chaosong/japan_chao/stationlist'), 'r');
format = '%s %f %f %f \n';
catcell = textscan(fid, format);
fclose(fid);
stnm = char(catcell{1,1});    % station name, string
stlo = catcell{1,2};    % station longitude
stla = catcell{1,3};    % station latitude
stel = catcell{1,4};    % station elevation

scatter(f.ax,stlo, stla,60,'y','filled','^','markeredgec','k');

text(f.ax,stlo, stla+0.03, stnm, 'fontsize',6);

% for test purpose, choose some stations
ind = [6;10;13;21;23;24;30;31];
stnmsel = stnm(ind,:);
stlosel = stlo(ind);
stlasel = stla(ind);
stelsel = stel(ind);

scatter(f.ax,stlosel, stlasel,60,'k','filled','^','markeredgec','k');

plot(f.ax,[135.62,136.22],[34.78,33.88],'color',...
        [0.6 0.6 0.6],'linew',2);


%% set the further boundary to regular events 
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

% events inside bound 1
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),bnd1(:,1),bnd1(:,2));
isinbnd1 = is | ion;
evtbnd1 = evtall(isinbnd1 == 1, :);

% events inside bound 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),bnd2(:,1),bnd2(:,2));
isinbnd2 = is | ion;
evtbnd2 = evtall(isinbnd2 == 1, :);

plot(f.ax,bnd2(:,1),bnd2(:,2),'color','b','linew',2);


%% create a cross-scetion profile, 
% csp1x = linspace(135.8,135.2,1000);
% csp1y = linspace(33.6,34.5,1000);
% csp1 = [csp1x' csp1y'];
% 
% csp2x = linspace(135.53,135.4,1000);
% csp2y = linspace(33.42,34.63,1000);
% csp2 = [csp2x' csp2y'];
% 

% cross-section region (rect)
csp = [135.5  33.4;
       135.38 34.6;
       135.2904 34.5910;
       135.4104 33.3910;
       135.5 33.4];

% the middle line of the cross-section
cspmidx = linspace(0.5*(csp(4,1)+csp(5,1)),0.5*(csp(2,1)+csp(3,1)),1000);
cspmidy = linspace(0.5*(csp(4,2)+csp(5,2)),0.5*(csp(2,2)+csp(3,2)),1000);
cspmid = [cspmidx' cspmidy'];

% cross-section region (rect)
csp2 = [135.95  33.6;
       135.83 34.8;
       135.7404 34.7910;
       135.8604 33.5910;
       135.96 33.6];   
   
% the middle line of the cross-section
csp2midx = linspace(0.5*(csp2(4,1)+csp2(5,1)),0.5*(csp2(2,1)+csp2(3,1)),1000);
csp2midy = linspace(0.5*(csp2(4,2)+csp2(5,2)),0.5*(csp2(2,2)+csp2(3,2)),1000);
csp2mid = [csp2midx' csp2midy'];

% assume tremors are on the slab interface, interpolate to get their depth
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
obara(:,8) = F(obara(:,6),obara(:,7));

patch(f.ax,csp(:,1),csp(:,2),'r','Facealpha',0.4);
plot(f.ax,cspmid(:,1),cspmid(:,2),'r','linew',1.5);
patch(f.ax,csp2(:,1),csp2(:,2),'r','Facealpha',0.4);
plot(f.ax,csp2mid(:,1),csp2mid(:,2),'r','linew',1.5);
text(f.ax,135.3,34.67,'A','FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,135.45,33.33,"A'",'FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,135.75,34.85,'B','FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
text(f.ax,135.85,33.55,"B'",'FontSize',10,'EdgeColor','k','Margin',1,'backgroundcolor','w');
title(f.ax,prefix);
f.ax.Box='on';

% save figure
print(f.fig,'-dpdf',strcat(figpath,'/evt_tmr_sta_bounds_csection.',prefix,'.dep',...
      num2str(depran),'.pdf'));   

  
%%  project all events & tremor onto it
% get events & tremors inside the profile, regardless of bound1 and 2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),csp(:,1),csp(:,2));
isincspe = is | ion;
evtcsp = evtall(isincspe == 1, :);

[is,ion] = inpolygon(obara(:,6),obara(:,7),csp(:,1),csp(:,2));
isincspt = is | ion;
tmrcsp = obara(isincspt == 1, :);

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
[dx,dy] = absloc2relaloc(evtcsp(:,8),evtcsp(:,9),lon0,lat0);
% rotate the coordinates to profile axis 
[evtcsp(:,13),evtcsp(:,14)] = coordinate_rot(dx,dy,rotang,0,0);

% convert tremor's locations to relative distance as well
[dx,dy] = absloc2relaloc(tmrcsp(:,6),tmrcsp(:,7),lon0,lat0);
[tmrcsp(:,9),tmrcsp(:,10)] = coordinate_rot(dx,dy,rotang,0,0);

% convert stations' locations to relative distance as well
[dx,dy] = absloc2relaloc(stlo,stla,lon0,lat0);
[stxr,styr] = coordinate_rot(dx,dy,rotang,0,0);


%%%% csp2
[is,ion] = inpolygon(evtall(:,8),evtall(:,9),csp2(:,1),csp2(:,2));
isincsp2e = is | ion;
evtcsp2 = evtall(isincsp2e == 1, :);

[is,ion] = inpolygon(obara(:,6),obara(:,7),csp2(:,1),csp2(:,2));
isincsp2t = is | ion;
tmrcsp2 = obara(isincsp2t == 1, :);

% define the origin as the mid lower right point
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


%% plot the cross-section profile
% x axis is along the profile, y axis is depth

f2.fig=figure;
widin = 10;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

% plot events inside csp
f2.ax(1) = subplot(221);
ax = f2.ax(1);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
axis(ax, 'equal');
plot(ax,cspmidr(:,5),cspmidr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp(:,10),tmrcsp(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp(evtcsp(:,11)<1, :);  
scatter(ax,tmp(:,14),tmp(:,10),4,[0.5 0.8 1],'filled','o');
% % tmp2 is event inside csp1, mag <1, and lfe
% tmptype = hypdet(isincspe == 1);
% tmp2 = evtcsp(evtcsp(:,11)<1 & (tmptype=='S' | tmptype=='s'), :);  
% scatter(ax,tmp2(:,14),tmp2(:,10),4,'g','filled','o');
scatter(ax,styr,zeros(length(styr),1),60,'y','filled','^','markeredgec','k');
text(ax,styr,5*ones(length(styr),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,80]);
xlim(ax,[-10,110]);
text(ax,60,70,"evts inside A-A', mag<1");
hold(ax,'off');

% plot events inside csp2
f2.ax(2) = subplot(222);
ax = f2.ax(2);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
axis(ax, 'equal');
plot(ax,csp2midr(:,5),csp2midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp2(:,10),tmrcsp2(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp2(evtcsp2(:,11)<1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,[0.5 0.8 1],'filled','o');
% % tmp2 is event inside csp1, mag <1, and lfe
% tmptype2 = hypdet(isincsp2e == 1);
% tmp2 = evtcsp2(evtcsp2(:,11)<1 & (tmptype2=='S' | tmptype2=='s'), :);
% scatter(ax,tmp2(:,14),tmp2(:,10),4,'g','filled','o');
scatter(ax,sty2r,zeros(length(sty2r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty2r,5*ones(length(sty2r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,80]);
xlim(ax,[-10,110]);
text(ax,60,70,"evts inside B-B', mag<1");
hold(ax,'off');

f2.ax(3) = subplot(223);
ax = f2.ax(3);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
axis(ax, 'equal');
plot(ax,cspmidr(:,5),cspmidr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp(:,10),tmrcsp(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp(evtcsp(:,11)>=1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,'b','filled','o');
% % tmp2 is event inside csp1, mag <1, and lfe
% tmp2 = evtcsp(evtcsp(:,11)>=1 & (tmptype=='S' | tmptype=='s'), :);
% scatter(ax,tmp2(:,14),tmp2(:,10),4,'g','filled','o');
scatter(ax,styr,zeros(length(styr),1),60,'y','filled','^','markeredgec','k');
text(ax,styr,5*ones(length(styr),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,80]);
xlim(ax,[-10,110]);
text(ax,60,70,"evts inside A-A', mag>=1");
hold(ax,'off');

f2.ax(4) = subplot(224);
ax = f2.ax(4);
hold(ax,'on');
grid(ax, 'on');
ax.Box = 'on';
axis(ax, 'equal');
plot(ax,csp2midr(:,5),csp2midr(:,3),'color','r','linew',0.5);     % slab along the profile
scatter(ax,tmrcsp2(:,10),tmrcsp2(:,8),4,[0.3 0.3 0.3],'filled','o');
tmp = evtcsp2(evtcsp2(:,11)>=1, :);
scatter(ax,tmp(:,14),tmp(:,10),4,'b','filled','o');
% % tmp2 is event inside csp1, mag <1, and lfe
% tmp2 = evtcsp2(evtcsp2(:,11)>=1 & (tmptype2=='S' | tmptype2=='s'), :);
% scatter(ax,tmp2(:,14),tmp2(:,10),4,'g','filled','o');
scatter(ax,sty2r,zeros(length(sty2r),1),60,'y','filled','^','markeredgec','k');
text(ax,sty2r,5*ones(length(sty2r),1),stnm,'fontsize',5);
set(ax, 'YDir','reverse');
set(ax, 'XDir','reverse');
ylabel(ax,{'Depth (km)'});
xlabel(ax,{'SE-NW along profile (km)'});
ylim(ax,[0,80]);
xlim(ax,[-10,110]);
text(ax,60,70,"evts inside B-B', mag>=1");
hold(ax,'off');

supertit(f2.ax,prefix);
% supertit(f2.ax,strcat(prefix,'--low/inferior precision hopocenter identification in green'));

% save figure
print(f2.fig,'-depsc2',strcat(figpath,'/evt_tmr_csection_view.',prefix,'.dep',...
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
view(ax,285,15);
zlim(ax,[20 80]);
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
zlim(ax,[20 80]);
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
zlim(ax,[20 80]);
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
zlim(ax,[20 80]);
hold off

supertit(f3.ax,prefix);

% save figure
savefig(f3.fig,strcat(figpath,'/evt_tmr_3dview.',prefix,'.dep',...
      num2str(depran),'.fig'));










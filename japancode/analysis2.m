% function analysis2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the advanced version main script to analyze the spatio-temporal tremor
% and regular events in the Japanese region (now in kii and shikoku) 
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
% First created date:   2020/02/24
% Last modified date:   2020/02/28
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

% how detailed to distinguish the magnitude difference, 
% class = 7;  % for a detailed magnitude difference
class = 3;  % only distinguish for all mag, mag<1, >=1 

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


%% spatial relationship between tremor & events inside & outside
f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
hold on
plot(rectbd(:,1),rectbd(:,2),'yo-','linew',1.5,'MarkerS',4);    
p1=scatter(obara(:,6),obara(:,7),6,[0.4 0.4 0.4],'filled','o');
p2=plot(tmrbd(:,1),tmrbd(:,2),'ko-','linew',1.5,'MarkerS',4);
p3=scatter(evtinall(:,8),evtinall(:,9),10,'b','filled','o');
p4=scatter(evtoutall(:,8),evtoutall(:,9),10,'r','filled','o');

legend([p1,p2,p3,p4],{'tremor','boundary','events inside','events outside'},'location','northwest');
axis equal
box on
grid on
% save figure
print(f.fig,'-dpdf',strcat(figpath,'/evt_tmr_with_bounds.',prefix,'.dep',num2str(depran),'.pdf'));


%% plot the occuring frequency(times)--magnitude relation of regular event 
f1.fig=figure;
widin = 5;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
f1.ax=pltevtmagfreq(evtinall(:,11),0.5);

% save figure
print(f1.fig,'-dpdf',strcat(figpath,'/evt_in_tmrbound.mag_freq.',prefix,'.dep',num2str(depran),'.pdf'));


%% plot the spatial distribution of the regular events with different magnitudes
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8.5;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end
% reposition
set(f2.ax(1), 'position', [ 0.08, 0.55, 0.4, 0.4], 'unit','normalized');
set(f2.ax(2), 'position', [ 0.56, 0.55, 0.4, 0.4], 'unit','normalized');
set(f2.ax(3), 'position', [ 0.08, 0.06, 0.4, 0.4], 'unit','normalized');
set(f2.ax(4), 'position', [ 0.56, 0.06, 0.4, 0.4], 'unit','normalized');

minmag = 1;
evtinm1 = evtinall(evtinall(:,11)>=minmag, :);

f2 = pltevtloc(f2,regflag,evtinm1,obara,slab);

% save figure
print(f2.fig,'-dpdf',strcat(figpath,'/evt_in_tmrbound.spatialloc.',prefix,'.dep',...
      num2str(depran),'.pdf'));


%% regular events occurence rate variation VS tremor detection rate, regardless of tremor ETS period
%%% for all events
yrall = unique(evtall(:,2));
% occurence rate, num per day of every year
nareaall = areaall/areaall;
rtayr = evtrateyr(yrall,evtall,obara,nareaall,class);

%%% for only events inside the tremor boundary
nareain = areain/areaall;
rtiyr = evtrateyr(yrall,evtinall,obara,nareain,class);

%%% for only events outside the tremor boundary
nareaout = (areaall-areain)/areaall;
rtoyr = evtrateyr(yrall,evtoutall,obara,nareaout,class);

%%% for all tremor detections
rttmryr = tmrrateyr(yrall,obara,nareain);

f3.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f3.ax(isub) = subplot(nrow,ncol,isub);
end

f3=pltevtrateyr(f3,yrall,rtayr,rtiyr,rtoyr,rttmryr,class);

% save figure
print(f3.fig,'-dpdf',strcat(figpath,'/evt_rate_yr.',num2str(class),'.',prefix,'.dep',...
      num2str(depran),'.pdf'));
  

%% regular event rate during tremor days & free days, yearly change.
datetmr = unique(obara(:,1));
yrall = unique(obara(:,2));

%%% rate for all events in rectangle during tremor days & free days
nareaall = areaall/areaall;
[rtatyr,rtafyr] = evtrateintmryr(yrall,evtall,obara,datetmr,nareaall,class);
        
%%% rate for only events inside tremor boundary during tremor days & free days
nareain = areain/areaall;
[rtityr,rtifyr] = evtrateintmryr(yrall,evtinall,obara,datetmr,nareain,class);

%%% rate for only events outside tremor boundary during tremor days & free days
nareaout = (areaall-areain)/areaall;
[rtotyr,rtofyr] = evtrateintmryr(yrall,evtoutall,obara,datetmr,nareaout,class);

f4.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 10.5;   % maximum height allowed is 11 inches
set(f4.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 2;
for isub = 1:nrow*ncol
    f4.ax(isub) = subplot(nrow,ncol,isub);
end
 
f4=pltevtratetmryr(f4,yrall,rtatyr,rtafyr,rtityr,rtifyr,rtotyr,rtofyr,class);

% save figure
print(f4.fig,'-dpdf',strcat(figpath,'/evt_rate_tmrdays_yr.',num2str(class),'.',prefix,'.dep',...
      num2str(depran),'.pdf'));


%% Think about ETS, not yearly change
%%% yearly change may not be representative so much, because ETS length, frequency, intensity might
%%% be different from year to year, so that the occurence rate averaged every year does not give
%%% enough information, so now i am using ETS period to confine the two types of time   
[stday, edday] = ETSestimate(regflag);

%%% rate for all events in rectangle during ETS or inter-ETS
nareaall = areaall/areaall;
[rtaetsm,rtanom] = evtrateets(stday,edday,evtall,obara,nareaall,class);
        
%%% rate for only events inside tremor boundary during ETS or inter-ETS
nareain = areain/areaall;
[rtietsm,rtinom] = evtrateets(stday,edday,evtinall,obara,nareain,class);

%%% rate for only events outside tremor boundary during ETS or inter-ETS
nareaout = (areaall-areain)/areaall;
[rtoetsm,rtonom] = evtrateets(stday,edday,evtoutall,obara,nareaout,class);

%%% rate for tremor detection during ETS
[rttmrets] = tmrrateets(stday,edday,obara,nareain);

f5.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 10.5;   % maximum height allowed is 11 inches
set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 3;
ncol = 2;
for isub = 1:nrow*ncol
    f5.ax(isub) = subplot(nrow,ncol,isub);
end
 
f5=pltevtrateets(f5,rtaetsm,rtanom,rtietsm,rtinom,rtoetsm,rtonom,rttmrets,class);

% save figure
print(f5.fig,'-dpdf',strcat(figpath,'/evt_rate_ets.',num2str(class),'.',prefix,'.dep',...
      num2str(depran),'.pdf'));


%% think about the event rate within a radius within a time range of a tremor detection
% imagine there is a tremor detection at a time at a location, what is the detection rate inside
% a horizontal radius within a time range relative to the tremor timing; in the same time range,
% what is the rate for the region outside that circle, but inside tremor region boundary; And what
% is the case for the region outside tremor region boundary

dttol = 2;      % max half time difference, in hr
dloctol = 0.1;  % max half location difference, in deg

% rate of events inside radius, outside radius, outside tremor boundary
% events inside radius, outside radius, outside tremor boundary
% tremors inside radius, outside radius, outside tremor boundary
[rteirm,rteorm,rteotm,evtir,evtor,evtot,tmrir,tmror,tmrot] = ...
         evtratehr(evtinall,evtoutall,obara,areaall,areain,dttol,dloctol,class);
     
% plot
f6.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f6.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 2;
for isub = 1:nrow*ncol
    f6.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f3.ax(1), 'position', [ 0.08, 0.56, 0.4, 0.4]);
set(f3.ax(2), 'position', [ 0.56, 0.56, 0.4, 0.4]);
set(f3.ax(3), 'position', [ 0.08, 0.08, 0.4, 0.4]);
set(f3.ax(4), 'position', [ 0.56, 0.08, 0.4, 0.4]);

minmag = -10;
% subplot 1, regular events vs. tremor during tremor active days
f6.ax(1) = pltevttmrloc2(f6.ax(1),regflag,slab,evtir,tmrir,minmag);

% subplot 2, regular events vs. tremor during tremor days, but inactive
f6.ax(2) = pltevttmrloc2(f6.ax(2),regflag,slab,evtor,tmror,minmag);

% subplot 3, regular events in tremor free days vs. all tremor
f6.ax(3) = pltevttmrloc2(f6.ax(3),regflag,slab,evtot,tmrot,minmag);

% subplot 4, regular events in tremor inactive days vs tremor in active days
hold(f6.ax(4),'on');
f6.ax(4).Box = 'on';
grid(f6.ax(4),'on');
corder=f6.ax(4).ColorOrder;
p1=scatter(f6.ax(4),1,rteirm(1,:),30,'ko','linew',1.5);
for i=2:size(rteirm,1)
    scatter(f6.ax(4),1,rteirm(i,:),25,corder(i,:),'o','linew',1.5);
end
p2=scatter(f6.ax(4),2,rteorm(1,:),30,'k^','linew',1.5);
for i=2:size(rteorm,1)
    scatter(f6.ax(4),2,rteorm(i,:),25,corder(i,:),'^','linew',1.5);
end
p3=scatter(f6.ax(4),3,rteotm(1,:),30,'kd','linew',1.5);
for i=2:size(rteotm,1)
    scatter(f6.ax(4),3,rteotm(i,:),25,corder(i,:),'d','linew',1.5);
end
ylabel(f6.ax(4),'Events per tremor detection per unit area');
xlim(f6.ax(4),[0 4]);
if class == 7
    p4 = scatter(f6.ax(4),1,rteirm(2,:),25,corder(2,:),'o','linew',1.5);
    p5 = scatter(f6.ax(4),1,rteirm(3,:),25,corder(3,:),'o','linew',1.5);
    p6 = scatter(f6.ax(4),1,rteirm(4,:),25,corder(4,:),'o','linew',1.5);
    p7 = scatter(f6.ax(4),1,rteirm(5,:),25,corder(5,:),'o','linew',1.5);
    p8 = scatter(f6.ax(4),1,rteirm(6,:),25,corder(6,:),'o','linew',1.5);
    p9 = scatter(f6.ax(4),1,rteirm(7,:),25,corder(7,:),'o','linew',1.5);
    ylim(f6.ax(4),[0 2]);
    legend(f6.ax(4),[p1,p2,p3,p4,p5,p6,p7,p8,p9],{'In radius, All mag','Out radius','Out boundary',...
        'Mag < 0','Mag = 0-1','Mag = 1-2','Mag = 2-3','Mag = 3-4', 'Mag >= 4'}, ...
        'location','northeast');
elseif class == 3
    p4 = scatter(f6.ax(4),1,rteirm(2,:),25,corder(2,:),'o','linew',1.5);
    p5 = scatter(f6.ax(4),1,rteirm(3,:),25,corder(3,:),'o','linew',1.5);
    ylim(f6.ax(4),[0 1.5]);
    legend(f6.ax(4),[p1,p2,p3,p4,p5],{'In radius, All mag','Out radius','Out boundary',...
        'Mag < 1', 'Mag >= 1'},'location','northeast');
end
hold(f6.ax(4),'off');

% save figure
print(f6.fig,'-dpdf',strcat(figpath,'/evt_tmr_hr_loc_rate.',num2str(class),'.',prefix,'.dep',...
      num2str(depran),'.pdf'));
     
            
















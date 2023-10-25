% function jma_vel_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to analyze the JMA2001 reference velocity model. Visualize
% and modify it to create an eligible input for taup_create, 
%
%
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/03/25
% Last modified date:   2020/03/25
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

%% load the velocity model
data = load('/home/data2/chaosong/shikoku_kii/jmavelmod/vjma2001');
vp = data(:,1);
vs = data(:,2);
dep = data(:,3);

% plot it      
f.fig=figure;
widin = 6;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

subplot(3,1,1);
ax = gca;
hold(ax,'on');
stairs(ax,vp,dep,'linew',0.5,'color','b');
stairs(ax,vs,dep,'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
text(0.83,0.91,'JMA2001','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel('Velocity (km/s)');
ylabel('Depth (km)');


subplot(3,1,2);
ax = gca;
hold(ax,'on');
stairs(ax,vp,dep,'linew',0.5,'color','b');
stairs(ax,vs,dep,'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
legend('Vp','Vs','location','best','fontsize',10);
ylim([0, 100]);
xlim([2, 8]);

subplot(3,1,3);
ax = gca;
hold(ax,'on');
stairs(ax,vp,dep,'linew',0.5,'color','b');
stairs(ax,vs,dep,'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
ylim([0, 20]);
xlim([2, 8]);


%% construct a new model 1, a densely sampled smooth and continuous model without any discontinuities 
fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ak135.tvel', 'r');
format = '%f %f %f %f \n';
dcell = textscan(fid,format,'headerlines',2); 
ak135(:,1) = dcell{1,1};    % dep
ak135(:,2) = dcell{1,2};    % vp
ak135(:,3) = dcell{1,3};    % vs

indak = find(ak135(:,1)>=max(dep),1,'first');
newmod = [dep vp vs; 
          ak135(indak:end,1) ak135(indak:end,2) ak135(indak:end,3)];

fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_test.tvel', 'w+');
fprintf(fid,'%s \n', 'jma2001-test - P');
fprintf(fid,'%s \n', 'jma2001-test - S');
fprintf(fid,'%10.4f %10.4f %10.4f\n',newmod');
fclose(fid);

% fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_dense.tvel', 'w+');
% fprintf(fid,'%s \n', 'jma2001-dense - P');
% fprintf(fid,'%s \n', 'jma2001-dense - S');
% fprintf(fid,'%10.4f %10.4f %10.4f\n',newmod');
% fclose(fid);

% plot it      
f.fig=figure;
widin = 6;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

subplot(2,1,1);
ax = gca;
hold(ax,'on');
plot(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
plot(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
text(0.75,0.91,'JMA2001_dense','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
ylim([0, 1000]);


subplot(2,1,2);
ax = gca;
hold(ax,'on');
plot(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
plot(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
legend('Vp','Vs','location','best','fontsize',10);
ylim([0, 100]);
xlim([2, 8]);


%% construct a new model 2, a sparsely sampled smooth and continuous model without any discontinuities 
fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ak135.tvel', 'r');
format = '%f %f %f %f \n';
dcell = textscan(fid,format,'headerlines',2); 
ak135(:,1) = dcell{1,1};    % dep
ak135(:,2) = dcell{1,2};    % vp
ak135(:,3) = dcell{1,3};    % vs

indak = find(ak135(:,1)>=max(dep),1,'first');
ind1 = find(dep>=4,1,'first');
ind2 = find(dep>=20,1,'first');
ind3 = find(dep>=30,1,'first');
ind4 = [];
for idep = 45: 1: 60    % sampling at 1 km interval
    tmp = find(dep>=idep,1,'first');
    ind4 = [ind4; tmp]; 
end
ind5 = find(dep>=200,1,'first');
ind6 = find(dep>=410,1,'first');
ind7 = [];
for idep = 600: 50: max(dep)   % sampling at 50 km interval, this depth should not affect much
    tmp = find(dep>=idep,1,'first');
    ind7 = [ind7; tmp]; 
end
indjma = [1; ind1; ind2; ind3; ind4; ind5; ind6; ind7];

newmod = [dep(indjma) vp(indjma) vs(indjma);
          ak135(indak:end,1) ak135(indak:end,2) ak135(indak:end,3)];

fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.tvel', 'w+');
fprintf(fid,'%s \n', 'jma2001-sparse - P');
fprintf(fid,'%s \n', 'jma2001-sparse - S');
fprintf(fid,'%10.4f %10.4f %10.4f\n',newmod');
fclose(fid);

% plot it      
f.fig=figure;
widin = 6;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

subplot(2,1,1);
ax = gca;
hold(ax,'on');
plot(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
plot(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
text(0.75,0.91,'JMA2001_sparse','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
ylim([0, 1000]);


subplot(2,1,2);
ax = gca;
hold(ax,'on');
plot(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
plot(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
legend('Vp','Vs','location','best','fontsize',10);
ylim([0, 100]);
xlim([2, 8]);


%% construct a new model 3, from Ukawa et al. 1983, with discontinuities 
fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ak135.tvel', 'r');
format = '%f %f %f %f \n';
dcell = textscan(fid,format,'headerlines',2); 
ak135(:,1) = dcell{1,1};    % dep
ak135(:,2) = dcell{1,2};    % vp
ak135(:,3) = dcell{1,3};    % vs

indak = find(ak135(:,1)>400,1,'first');

% The following is from Ukawa etal 1983, table 2 on page 10
ukawa83 = [0   5.50 3.25;   
           10  5.98 3.49;
           20  6.51 3.74;
           32  7.20 4.07;
           32  7.80 4.41;
           50  7.85 4.43;
           100 8.00 4.50;
           150 8.12 4.53;
           200 8.26 4.60;
           250 8.42 4.68;
           300 8.58 4.76;
           350 8.75 4.85];
%            400 8.92 4.94];    % neglect the 400 km, otherwise it would causes a vel decrease there

newmod = [ukawa83;
          ak135(indak:end,1) ak135(indak:end,2) ak135(indak:end,3)];

fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ukawa83.tvel', 'w+');
fprintf(fid,'%s \n', 'ukawa83 - P');
fprintf(fid,'%s \n', 'ukawa83 - S');
fprintf(fid,'%10.4f %10.4f %10.4f\n',newmod');
fclose(fid);

% plot it      
f.fig=figure;
widin = 6;  % maximum width allowed is 8.5 inches
htin = 10;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/30 widin*res htin*res]);

subplot(2,1,1);
ax = gca;
hold(ax,'on');
plot(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
plot(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
text(0.75,0.91,'Ukawa83','FontSize',12,'unit','normalized','EdgeColor','k','Margin',2);
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
ylim([0, 1000]);


subplot(2,1,2);
ax = gca;
hold(ax,'on');
plot(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
plot(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
xlabel('Velocity (km/s)');
ylabel('Depth (km)');
legend('Vp','Vs','location','best','fontsize',10);
ylim([0, 100]);
xlim([2, 8]);


%% use taupcreate to make a new taup-style model for further use
% In default, taupcreate in mattaup creat the model in the path same as the mfile, the .taup
% velocity model can ONLY been called by absolute path;
% This is different from the linux version taup_create, which can be called as long as the model is
% placed into the StdModel folder. And models from the two seem UNcompatible with each other. 
% '.tvel' file can be read by linux and matlab, but you need to create '.taup' separately in
% different platforms

%%%% for jma2001 model
intvel1 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.tvel';
outtaup1 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.taup';
taupcreate(intvel1, outtaup1);    

tic
tt=tauptime('mod',outtaup1,'dep',50,'ph','s,S','dist',100);
toc

%%%% for ukawa83 model
intvel2 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ukawa83.tvel';
outtaup2 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ukawa83.taup';
taupcreate(intvel2, outtaup2);    

tic
tt=tauptime('mod',outtaup2,'dep',50,'ph','s,S','dist',100);
toc

%% test
intvel1 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_test.tvel';
outtaup1 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_test.taup';
taupcreate(intvel1, outtaup1);    

tic
tt=tauptime('mod',outtaup1,'dep',50,'ph','s,S','dist',100);
toc


















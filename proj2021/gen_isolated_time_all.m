% function gen_isolated_time_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to obtain the corresponding time offset from the 
% representative isolated detections from the rough identification algorithm.
%
% We aim to find the location relationship between these 'uncommon' detections
% relative to the HF locations or typical LF detections in which lower-freq
% energy is there, but the amplitude is smaller than the that of the high-freq
% energy. This would give us some hints on attenuation effect
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/06/14
% Last modified date:   2021/06/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

% set path
workpath = getenv('MHOME');
pgcpath = strcat(workpath, '/Seisbasics/hypoinverse/forproj21');

%% PGC trio
% relative arrival time to main station for each fam
contpgcoff = [
              86 20;    % fam 002
              61 29;  % fam 243 
              37 -4;  % fam 240
              ];
          
%%% load the template filtering effect result
workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

fam='002';     % family number
lolf = 0.5;
hilf = 1.25;
lohf = 1.25;
hihf = 6.5;
winsechf = 4;
winseclf = 16;
loopoffmax = 4;
xcmaxAVEnmin = 0.5;
sps = 20;

PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
    '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(loopoffmax),'.ccm', ...
    num2str(xcmaxAVEnmin),'.','80sps');
fname = strcat(rstpath, '/MAPS/tempfeff_PGC_enc','_',PREFIX);
filtcor = load(fname);
spsratio = sps/80;
filtlf002 = filtcor(:,2)*spsratio;   % lf template shift due to filterig effect, sign is the same
filtlf002 = filtlf002';

fam='243';     % family number
lolf = 0.5;
hilf = 1.25;
lohf = 1.25;
hihf = 6.5;
winsechf = 4;
winseclf = 16;
loopoffmax = 4;
xcmaxAVEnmin = 0.45;
sps = 20;

PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
    '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(loopoffmax),'.ccm', ...
    num2str(xcmaxAVEnmin),'.','80sps');
fname = strcat(rstpath, '/MAPS/tempfeff_PGC_enc','_',PREFIX);
filtcor = load(fname);
spsratio = sps/80;
filtlf243 = filtcor(:,2)*spsratio;   % lf template shift due to filterig effect, sign is the same
filtlf243 = filtlf243';

fam='240';     % family number
lolf = 0.5;
hilf = 1.25;
lohf = 1.25;
hihf = 6.5;
winsechf = 4;
winseclf = 16;
loopoffmax = 4;
xcmaxAVEnmin = 0.45;
sps = 20;

PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
    '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(loopoffmax),'.ccm', ...
    num2str(xcmaxAVEnmin),'.','80sps');
fname = strcat(rstpath, '/MAPS/tempfeff_PGC_enc','_',PREFIX);
filtcor = load(fname);
spsratio = sps/80;
filtlf240 = filtcor(:,2)*spsratio;   % lf template shift due to filterig effect, sign is the same
filtlf240 = filtlf240';


arr2 = [];
% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution 
off002 = [
          14 0;
          0 -8;
          ];
off002 = off002 - filtlf002;

%offset in broader-band, 0.5-6.5 hz, since the filtering effect is obtained in the same passband,
%the alignment is f-e-free
off002bb = [
          13 0;
          0 -6;
          ];

tmp = [];
for i = 1: size(off002,1)
    tmp(i,1:2) = contpgcoff(1,:);
end
tmp(:, 3:4) = off002;
arr2 = [arr2; tmp];
    
% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution 
off243 = [
          -11 6;
          1 5;
          1 4;
          -8 1;
          -12 4;
          1 5;
          3 2;
          2 1;  % not a perfect one, but close in space and time to the one above
          ];
off243 = off243 - filtlf243;

%offset in broader-band, 0.5-6.5 hz, since the filtering effect is obtained in the same passband,
%the alignment is f-e-free
off243bb = [
          -12 5;
          -1 3;
          -1 4;
          -8 2;
          -9 5;
          -1 3;
          2 0;
          1 0;  % not a perfect one, but close in space and time to the one above
          ];

tmp = [];
for i = 1: size(off243,1)
    tmp(i,1:2) = contpgcoff(2,:);
end
tmp(:, 3:4) = off243;
arr2 = [arr2; tmp];
      
% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution 
off240 = [
          -5 0;
          -17 -19;
          ];
off240 = off240 - filtlf240;

%offset in broader-band, 0.5-6.5 hz, since the filtering effect is obtained in the same passband,
%the alignment is f-e-free
off240bb = [
          -5 2;
          -15 -19;
          ];

tmp = [];
for i = 1: size(off240,1)
    tmp(i,1:2) = contpgcoff(3,:);
end
tmp(:, 3:4) = off240;
arr2 = [arr2; tmp];
      
      
fid2 = fopen(strcat(pgcpath, '/offset_all_isolatedlf'),'w+');
fprintf(fid2,'%d %d %f %f \n',...
        arr2');
fclose(fid2);








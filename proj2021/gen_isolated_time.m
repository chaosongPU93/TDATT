% function gen_isolated_time
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
              86 20;
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
filtlf = filtcor(:,2)*spsratio;   % lf template shift due to filterig effect, sign is the same
filtlf = filtlf';

% +2/-2 samples in 40 sps, equivalent to +1/-1 samples in 20 sps, ie, LF resolution 
off2 = [11 0;
        11 1;
        12 -1;
        13 -3;
        -2 -8;
        -1 -8;
        -3 -9;
        0 -4
        ];   
off2 = off2 - filtlf;

arr2 = [];
for i = 1: size(off2,1)
    tmp(i,1:2) = contpgcoff;
end
tmp(:, 3:4) = off2;
arr2 = [arr2; tmp];
    

fid2 = fopen(strcat(pgcpath, '/offset_002_isolatedlf'),'w+');
fprintf(fid2,'%d %d %f %f \n',...
        arr2');
fclose(fid2);








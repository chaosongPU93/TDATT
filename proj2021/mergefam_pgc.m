%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to merge the detections from all families IF you want 
% to set a distance cutoff first, then merge fams to remove duplicates base
% on time and CC coef, so run this AFTER hypoinverse.
% The input is already been removed off the double counting. Since for each
% fam, it wiped put the same data, so the detections should still have a lot of
% duplicates, preserve the ones within the one duration that have the highest 
% CC value. In the mean time, detections in each family have a different
% reference (0,0) point, so the merging can also convert all to the same
% frame (similar to changing centroid).
%
% this is the version for merging fams using PGC trio, the idea behind is the 
% same as in 'merge_fam_afthypo_allfam' and 'merge_fam_afthypo_requestfam',
% depending on if the input data has been through distance cutoff
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2021/06/24
% Last modified date:   2021/06/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forproj21');

nfampool = ['002';
            '243';
            '240';
            '253';
            '036';
            '251';
            ];

nfam = size(nfampool,1);
disp(nfam); 

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
    
stas=['PGC  '
      'SSIB '
      'SILB '];     % determine the trio and order, here the 1st sta is PGC
        
winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.45;


%% load detections with or without a distance cutoff
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));

% dcutflag = 1;
dcutflag = 0;

if dcutflag
    distmaxhf = 10;
    distmaxlf = 10;
    hffname = strcat(rstpath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.pgcfam.pj21.',num2str(distmaxlf),'kmdcut.',SUFFIXlf);
else
    hffname = strcat(rstpath, '/evtloc.pgcfam.pj21.nodcut.',SUFFIXhf);
    lffname = strcat(rstpath, '/evtloc.pgcfam.pj21.nodcut.',SUFFIXlf);
end

hfold = load(hffname);
lfold = load(lffname);
% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep + previous 50 cols
%%% UPDATED at 2021/06/23
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+4;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)


%% sort according to time, remove double counting
disp('Start removing double counts ...');

%%% for hf
%%% sort according to date and main arrival time, could add offset
hfsort = sortrows(hfold, [14,16]);
dateall = unique(hfsort(:,14));
hfnew = [];
for id = 1: length(dateall)
    hfday = hfsort(hfsort(:,14) == dateall(id), :);
    dtminhf = 0.5;      % min time during which only 1 detection is retained
%     colnum = [16 17 55 56 57 58];
    colnum = [16 17];
    indexthf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    wlensechf = 4;
%     [hfdaynew,indsave,inddis] = RemoveDoubleCounting3(hfday,dtminhf,indexthf,colnum,wlensechf);
    [hfdaynew,indsave,inddis] = RemoveDoubleCounting(hfday,dtminhf,indexthf,colnum);
    hfnew = [hfnew; hfdaynew];    % hfdaynew is the file of all fam for one day  
    
end

%%% for lf
%%% sort according to date and main arrival time, could add offset
lfsort = sortrows(lfold, [14,16]);
dateall = unique(lfsort(:,14));
lfnew = [];
for id = 1: length(dateall)
    lfday = lfsort(lfsort(:,14) == dateall(id), :);
    dtminlf = 1.1;      % min time during which only 1 detection is retained
%     colnum = [16 17 55 56 57 58];
    colnum = [16 17];
    indextlf = length(nfampool);    % the extreme case is that every fam contains a duplicate
    wlenseclf = 16;
%     [lfdaynew,indsave,inddis] = RemoveDoubleCounting3(lfday,dtminlf,indextlf,colnum,wlenseclf);    
    [lfdaynew,indsave,inddis] = RemoveDoubleCounting(lfday,dtminlf,indextlf,colnum);    
    lfnew = [lfnew; lfdaynew];    % hfdaynew is the file of all fam for one day  
    
end

% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep + previous 50 cols
%%% UPDATED at 2021/06/23
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+4;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

if dcutflag 
    fidhf = fopen(strcat(rstpath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',SUFFIXhf),'w+');
    fidlf = fopen(strcat(rstpath, '/evtloc.pgcfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',SUFFIXlf),'w+');
else
    fidhf = fopen(strcat(rstpath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXhf),'w+');
    fidlf = fopen(strcat(rstpath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXlf),'w+');
end

fprintf(fidhf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        hfnew');
fclose(fidhf);

fprintf(fidlf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        lfnew');
fclose(fidlf);








function [trange, hfnew, lfnew] = catalog_of_bursts_in_region(FLAG,TYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT recommended any more
% Use 'catalog_of_region_of_interest' ---> 'catalog_of_tremor_bursts'
%
% --Problems come up when first running the 'catalog_of_tremor_bursts' and
%   and then the 'catalog_of_region_of_interest', which identify the tremor 
%   bursts of the entire catalog and generate a smaller catalog of bursts
%   only first, then restrain the detections in bursts to be inside some 
%   region and adjust the new time range of the burst.
% --Ideally we want the time ranges when tremor was active at different 
%   regions to be temporally exclusive, so that when we get the spectra of
%   these time windows, we are actually focusing on the region of interest.
% --But it seems not the case, the time ranges at diff regions are highly
%   overlapping, no matter if you use the max/min time of the in-bound 
%   detections to define the new time range or the prctile of the times
%   to avoid extremes.
% --3 possibilities might be there: Less likely is that the codes are doing
%   what i want, ie. getting the times correctly, NEEDs check!
%   Another one is that during most bursts of time of the entire catalog,
%   tremor was active over a larger region than your request, so it can't be 
%   temporally divided properly even you force a smaller region.
%   The final one may be is the order of running the codes is not ideal, ie.
%   if you specify the region and catalog of region first, then use the time
%   clustering algorithm to get the tremor bursts particularly inside the
%   region, the problem could be solved.
% --This is the purpose of this code
%
% --So the order would be opposite to the previous, it will
%   narrow down the region first, then get the bursts. This might lead to a
%   increase of inter-event time threshold and drop of min. number of detections
%   in the burst.
%   'catalog_of_region_of_interest' ---> 'catalog_of_tremor_bursts'
% --Note that the revised version only outputs the refined tremor catalog in
%   in region of interest from the raw, but NOT the time range anymore
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/10/26
% Last modified date:   2021/10/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
% clear
close all

%% default value for easy debugging
defval('FLAG','PGC');
defval('TYPE','HF');

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
hypopath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forproj21');

% FLAG = 'PGC';
% FLAG = 'TWKB';

disp(FLAG);

if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  temppath = strcat(datapath, '/templates/LZBtrio/');
  rstpath = strcat(datapath, '/LZBtrio');
  figpath = strcat(workpath, '/project2021/TWKBtrio');
  
  nfampool = [
              '002';
              '043';
              '141';
              '047';
              '010';
              '144';
              '099';
              '068';
              '125';
              '147';
              '017';
              '006';
              '001';
              ];
  nfam = size(nfampool,1);
  
  stas=['TWKB '
    'LZB  '
    'MGCB '];     % determine the trio and order
  nsta=size(stas,1);         %  number of stations
  
  %%% load detections with or without a distance cutoff
  winlensechf = 4;
  winlensec=16;
  loopoffmax=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
  xcmaxAVEnmin=0.35; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
  SUFFIXhf = strcat('up.hf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
    '.',num2str(xcmaxAVEnmin));
  SUFFIXlf = strcat('lf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
    '.',num2str(xcmaxAVEnmin));
%     hffname = strcat(hypopath, '/evtloc.lzbfam.pj21.nodcutnodou.',SUFFIXhf);
%     lffname = strcat(hypopath, '/evtloc.lzbfam.pj21.nodcutnodou.',SUFFIXlf);
  
  distmaxhf = 8;
  distmaxlf = 12;
%   if strcmp(CATA, 'raw')	% option 1: the raw merged catalog
    hffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',SUFFIXhf);
    lffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',SUFFIXlf);
%   elseif strcmp(CATA, 'burst')  % option 2: the refined catalog based on tremor bursts
%     hffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',...
%       TYPE,'.',SUFFIXhf);
%     lffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',...
%       TYPE,'.',SUFFIXlf);
%   end
  
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
  temppath = strcat(datapath, '/templates/PGCtrio/');
  rstpath = strcat(datapath, '/PGCtrio');
  figpath = strcat(workpath, '/project2021/PGCtrio');
  
 %     nfampool = [
 %                 '002'
 %                 '243';
 %                 '253';
 %                 '036';
 %                 '034';
 %                 '061';
 %                 '023';
 %                 '251';
 %                 '240';
 %                 '255';
 %                 '012';
 %                 '065';
 %                 ];
  
  nfampool = [
                '002';
                '243';
                '240';  % the most recent catalog of PGC is from fam 002, 243 and 240
%                 '253';
%                 '036';
%                 '251';
              ];
  nfam = size(nfampool,1);
  
  stas=['PGC  '
    'SSIB '
    'SILB '];
  nsta=size(stas,1);         %  number of stations
  
  %%% load the merged catalog of hf and lf detections
  winlensechf = 4;
  winlensec = 16;
  loopoffmax = 4;
  xcmaxAVEnmin = 0.45;
  SUFFIXhf = strcat('up.hf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',...
    num2str(loopoffmax),'.',num2str(xcmaxAVEnmin));
  SUFFIXlf = strcat('lf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
    '.',num2str(xcmaxAVEnmin));
  %     hffname = strcat(hypopath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXhf);
  %     lffname = strcat(hypopath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXlf);
  
  distmaxhf = 10;
  distmaxlf = 10;
%   if strcmp(CATA, 'raw')	% option 1: the raw merged catalog
    hffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',SUFFIXhf);
    lffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',SUFFIXlf);
%   elseif strcmp(CATA, 'burst')  % option 2: the refined catalog based on tremor bursts
%     hffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcutnodou.burst',...
%       TYPE,'.',SUFFIXhf);
%     lffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxlf),'kmdcutnodou.burst',...
%       TYPE,'.',SUFFIXlf);
%   end

end

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
  'LZB'];
POLSTA=['SSIB '           % polaris station names
  'SILB '
  'KLNB '
  'MGCB '
  'TWKB '];

hftime = load(hffname);
lftime = load(lffname);
% NOTE: if using the burst catalog, there are 59 cols where the last col is the occurrence time in
% days w.r.t. the starting day of each ETS, 2003060, 2004194, 2005254
%
% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep  (8 + previous 50 cols)
%%% UPDATED at 2021/06/23
%%%   n=n+8;
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
  
% sort according to day, sec of the window center
daycol = 14;
seccol = 15;
hftime = sortrows(hftime, [daycol, seccol]);
lftime = sortrows(lftime, [daycol, seccol]);

ncol = size(hftime,2);

% for 2005 only, rotate to the strike N30W and dip N60E direction, 2003 and 2004 have the main front
% migration direction to N
rotang = 30;  % counter-clockwise from y to strike, 2005

% keyboard

%% plot the cumulative density map to see the active regions
% 2003
hfplt = hftime(hftime(:,daycol) < 2004*1000, :);
lfplt = lftime(lftime(:,daycol) < 2004*1000, :);
xran = [-25 30];
yran = [-20 30];
[f1] = plt_cumulative_density(hfplt,lfplt,xran,yran);
text(f1.ax(1), 0.85, 0.93, '2003','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);

% 2004
hfplt = hftime(hftime(:,daycol) > 2004*1000 & hftime(:,daycol) < 2005*1000, :);
lfplt = lftime(lftime(:,daycol) > 2004*1000 & lftime(:,daycol) < 2005*1000, :);
xran = [-25 30];
yran = [-20 30];
[f2] = plt_cumulative_density(hfplt,lfplt,xran,yran);
text(f2.ax(1), 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);

% 2005
hfplt = hftime(hftime(:,daycol) > 2005*1000, :);
lfplt = lftime(lftime(:,daycol) > 2005*1000, :);
xran = [-25 30];
yran = [-20 30];
[f3] = plt_cumulative_density(hfplt,lfplt,xran,yran);
text(f3.ax(1), 0.85, 0.93, '2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);
      
% All
hfplt = hftime;
lfplt = lftime;
xran = [-25 30];
yran = [-20 30];
[f4] = plt_cumulative_density(hfplt,lfplt,xran,yran);
text(f4.ax(1), 0.85, 0.93, 'All','FontSize',12,'unit','normalized','horizontalalignment','center',...
     'EdgeColor','k','Margin',2);


%% define a region of interest
%%% this should come from the most activated region, or the tele-seismic receiver function (RF) study
% a region around fam 002 in Allan's 2003 paper 
bnd = [-10 -5;
       5 -5;
       5 5;
       -10 5;
       -10 -5;
      ];
hold(f4.ax(1),'on');
plot(f4.ax(1),bnd(:,1),bnd(:,2),'k-','linew',2);
hold(f4.ax(1),'off');
hold(f4.ax(2),'on');
plot(f4.ax(2),bnd(:,1),bnd(:,2),'k-','linew',2);
hold(f4.ax(2),'off');

% handpicked polygon 1 that is the SE part of the black box 
bndse = [3.0000e-01  -1.1000e+00
         1.1000e+00  -7.0000e-01
         2.1000e+00   3.0000e-01
         3.3000e+00   1.1000e+00
         4.5000e+00   1.9000e+00
         4.7000e+00  -7.0000e-01
         4.3000e+00  -1.7000e+00
         3.1000e+00  -2.9000e+00
         1.1000e+00  -3.5000e+00
         -7.0000e-01  -4.1000e+00
         -1.1000e+00  -3.1000e+00
         -5.0000e-01  -1.7000e+00
         3.0000e-01  -1.1000e+00
         ];
hold(f4.ax(1),'on');
plot(f4.ax(1),bndse(:,1),bndse(:,2),'k-','linew',2);
hold(f4.ax(1),'off');
hold(f4.ax(2),'on');
plot(f4.ax(2),bndse(:,1),bndse(:,2),'k-','linew',2);
hold(f4.ax(2),'off');

% handpicked polygon 2 that is the center part of the black box 
bndcnt = [-7.0000e-01   1.9000e+00
          5.0000e-01   2.1000e+00
          1.5000e+00   2.5000e+00
          2.5000e+00   2.3000e+00
          2.5000e+00   1.3000e+00
          1.7000e+00   3.0000e-01
          1.3000e+00  -3.0000e-01
          9.0000e-01  -9.0000e-01
          3.0000e-01  -1.3000e+00
          -9.0000e-01  -1.5000e+00
          -2.3000e+00  -2.1000e+00
          -2.7000e+00  -1.1000e+00
          -2.3000e+00   1.0000e-01
          -1.7000e+00   9.0000e-01
          -7.0000e-01   1.9000e+00
          ];
hold(f4.ax(1),'on');
plot(f4.ax(1),bndcnt(:,1),bndcnt(:,2),'k-','linew',2);
hold(f4.ax(1),'off');
hold(f4.ax(2),'on');
plot(f4.ax(2),bndcnt(:,1),bndcnt(:,2),'k-','linew',2);
hold(f4.ax(2),'off');

% handpicked polygon 3 that is the NW part of the black box , close to fam 047
bndnw = [-6.9000e+00  -9.0000e-01
         -6.5000e+00   1.0000e-01
         -6.1000e+00   1.5000e+00
         -4.7000e+00   2.5000e+00
         -3.7000e+00   3.3000e+00
         -2.5000e+00   4.1000e+00
         -7.0000e-01   4.5000e+00
         7.0000e-01   4.3000e+00
         1.0000e-01   3.1000e+00
         -7.0000e-01   1.9000e+00
         -1.5000e+00   1.1000e+00
         -2.5000e+00   5.0000e-01
         -3.7000e+00  -7.0000e-01
         -4.7000e+00  -1.7000e+00
         -6.1000e+00  -2.3000e+00
         -6.9000e+00  -9.0000e-01
         ];
hold(f4.ax(1),'on');
plot(f4.ax(1),bndnw(:,1),bndnw(:,2),'k-','linew',2);
hold(f4.ax(1),'off');
hold(f4.ax(2),'on');
plot(f4.ax(2),bndnw(:,1),bndnw(:,2),'k-','linew',2);
hold(f4.ax(2),'off');
       
% a rectangle N of the main box, close to fam 243
bndn = [-3 16;
       	-3 7;
        12 7;
        12 16;
        -3 16;
        ];
      
hold(f4.ax(1),'on');
plot(f4.ax(1),bndn(:,1),bndn(:,2),'k-','linew',2);
hold(f4.ax(1),'off');
hold(f4.ax(2),'on');
plot(f4.ax(2),bndn(:,1),bndn(:,2),'k-','linew',2);
hold(f4.ax(2),'off');

% a rectangle NE of the main box, close to fam 240
bndne = [12 18;
       	 12 8;
         18 8;
         18 18;
         12 18;
         ];
      
hold(f4.ax(1),'on');
plot(f4.ax(1),bndne(:,1),bndne(:,2),'k-','linew',2);
hold(f4.ax(1),'off');
hold(f4.ax(2),'on');
plot(f4.ax(2),bndne(:,1),bndne(:,2),'k-','linew',2);
hold(f4.ax(2),'off');

%% output the catalog of detections inside the region and time ranges
bnduse = bnd;

% step 1, events further inside the bound
[iin,ion] = inpolygon(hftime(:,1),hftime(:,2),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
hfbnd = hftime(isinbnd == 1, :);

[iin,ion] = inpolygon(lftime(:,1),lftime(:,2),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
lfbnd = lftime(isinbnd == 1, :);
    
% keyboard

% 
% if isequal(bnduse,bnd)
%   bndflag = 'bnd';
% elseif isequal(bnduse,bndse)
%   bndflag = 'bndse';
% elseif isequal(bnduse,bndcnt)
%   bndflag = 'bndcnt';
% elseif isequal(bnduse,bndnw)
%   bndflag = 'bndnw';
% elseif isequal(bnduse,bndn)
%   bndflag = 'bndn';
% elseif isequal(bnduse,bndne)
%   bndflag = 'bndne';
% end
% 
% if strcmp(FLAG, 'TWKB')     % use TWKB catalog
%   hffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxhf),'dcutnodb.',...
%     TYPE,bndflag,'.',SUFFIXhf);
%   lffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxlf),'dcutnodb.',...
%     TYPE,bndflag,'.',SUFFIXlf);  
% elseif strcmp(FLAG, 'PGC')     % use PGC catalog
%   hffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxhf),'dcutnodb.',...
%     TYPE,bndflag,'.',SUFFIXhf);
%   lffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxlf),'dcutnodb.',...
%     TYPE,bndflag,'.',SUFFIXlf);
% end
% 
% fidhf = fopen(hffname,'w+');
% fprintf(fidhf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f %.5f \n',...
%         hfnew');
% fclose(fidhf);
% 
% fidlf = fopen(lffname,'w+');
% fprintf(fidlf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f %.5f \n',...
%         lfnew');
% fclose(fidlf);
% 
% fid = fopen(strcat(hypopath, '/tremor_burst_ranges_',FLAG,'.',TYPE,'_',bndflag),'w+');
% fprintf(fid,'%d %9.2f %9.2f \n',ntrange');
% fclose(fid);
% 
% keyboard


%% separation in time between itself and its preceding
% obtain the relative time in days
ih03 = find(hfbnd(:,daycol) < 2004*1000);
hfbnd(ih03,ncol+1) = (hfbnd(ih03,daycol)-2003060)+hfbnd(ih03,seccol)./(3600.*24);
ih04 = find(hfbnd(:,daycol) < 2005*1000 & hfbnd(:,daycol) > 2004*1000);
hfbnd(ih04,ncol+1) = (hfbnd(ih04,daycol)-2004194)+hfbnd(ih04,seccol)./(3600.*24);
ih05 = find(hfbnd(:,daycol) > 2005*1000);
hfbnd(ih05,ncol+1) = (hfbnd(ih05,daycol)-2005254)+hfbnd(ih05,seccol)./(3600.*24);
% % rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% % to down-dip and strike
% [hfbnd(ih05,1),hfbnd(ih05,2)] = coordinate_rot(hfbnd(ih05,1),hfbnd(ih05,2),rotang,0,0);

il03 = find(lfbnd(:,daycol) < 2004*1000);
lfbnd(il03,ncol+1) = (lfbnd(il03,daycol)-2003060)+lfbnd(il03,seccol)./(3600.*24);
il04 = find(lfbnd(:,daycol) < 2005*1000 & lfbnd(:,daycol) > 2004*1000);
lfbnd(il04,ncol+1) = (lfbnd(il04,daycol)-2004194)+lfbnd(il04,seccol)./(3600.*24);
il05 = find(lfbnd(:,daycol) > 2005*1000);
lfbnd(il05,ncol+1) = (lfbnd(il05,daycol)-2005254)+lfbnd(il05,seccol)./(3600.*24);
% % rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% % to down-dip and strike
% [lfbnd(il05,1),lfbnd(il05,2)] = coordinate_rot(lfbnd(il05,1),lfbnd(il05,2),rotang,0,0);

% for the detections, obtain the separation in time between itself and its preceding
% detection for the catalog of each ETS.
% 1st col: occurence time
% 2nd col: inter-event time to its preceding detection
hfinter03 = interevt_time(hfbnd(ih03,end));
hfinter04 = interevt_time(hfbnd(ih04,end));
hfinter05 = interevt_time(hfbnd(ih05,end));
hfinter = [hfinter03; hfinter04; hfinter05];

% keyboard

%% plot the inter-detection time
%%% HF catalog
% set a threshold of inter-detection time
ttol1 = 1e-3;
% ttol2 = 1e-3;
ttol2 = 8*mad(hfinter03(:,2), 1);
ymax = 2e-3;
[f] = plt_interevt_time_diffets(hfinter03,hfinter04,hfinter05,ymax,ttol1,ttol2);
text(f.ax(1),0.04,0.9,'HF','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% print(f.fig,'-dpdf',strcat(figpath,'/inter-detection_time_hf_',FLAG,'.pdf'));  

% keyboard

%% group the tremor bursts 
%%% group the indice of the tremor bursts accroding to the threshold on the inter-detection time and
%%% number of detections in the burst for both HF and LF

if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  % set a threshold of inter-detection time
  ttolhf = 1e-3;
%   ttolhf = 5e-4;
  % set a threshold of minimum number of detections in the burst period
%   ntolhf = 15;
%   ntollf = 15;
  ntolhf = 30;
  ntollf = 15;
  [bursthf03, nburstlf03] = group_tremor_burst(hfinter03,lfbnd(il03,end),ttolhf,ntolhf,0);
  [bursthf04, nburstlf04] = group_tremor_burst(hfinter04,lfbnd(il04,end),ttolhf,ntolhf,0);
  [bursthf05, nburstlf05] = group_tremor_burst(hfinter05,lfbnd(il05,end),ttolhf,ntolhf,0);
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
%   ttolhf = 1e-3;  % param used until 2021/10/24
%   ntolhf = 30;  % param used until 2021/10/24
%   ttolhf = 2e-3;  % param used from 2021/10/25
%   ntolhf = 20;  % param used from 2021/10/25
  ttolhf = ttol2;  % param used from 2021/10/27, using MAD statistics 
  ntolhf = 15;  % param used from 2021/10/27
  ntollf = 10;
  [bursthf03, nburstlf03] = group_tremor_burst(hfinter03,lfbnd(il03,end),ttolhf,ntolhf,0);
  [bursthf04, nburstlf04] = group_tremor_burst(hfinter04,lfbnd(il04,end),ttolhf,ntolhf,0);
  [bursthf05, nburstlf05] = group_tremor_burst(hfinter05,lfbnd(il05,end),ttolhf,ntolhf,0);
end

% recover the index of the detections to occurence times with the same format as trange
if ~isempty(bursthf03)
  [thf03,tranhf03,nhf03] = burst_range(bursthf03,hfinter03,2003060);
  perchf03 = nhf03/length(ih03)*100;
  nlf03 = sum(nburstlf03);
  perclf03 = nlf03/length(il03)*100;
else
  thf03 = []; tranhf03 = []; nhf03 = 0; perchf03 = 0; nlf03 = 0;
end

if ~isempty(bursthf04)
  [thf04,tranhf04,nhf04] = burst_range(bursthf04,hfinter04,2004194);
  perchf04 = nhf04/length(ih04)*100;
  nlf04 = sum(nburstlf04);
  perclf04 = nlf04/length(il04)*100;
else
  thf04 = []; tranhf04 = []; nhf04 = 0; perchf04 = 0; nlf04 = 0;
end

if ~isempty(bursthf05)
  [thf05,tranhf05,nhf05] = burst_range(bursthf05,hfinter05,2005254);
  perchf05 = nhf05/length(ih05)*100;
  nlf05 = sum(nburstlf05);
  perclf05 = nlf05/length(il05)*100;
else
  thf05 = []; tranhf05 = []; nhf05 = 0; perchf05 = 0; nlf05 = 0;
end

perchf = (nhf03+nhf04+nhf05)/(length(ih03)+length(ih04)+length(ih05))*100;

% keyboard


%% plot the selected tremor bursts
if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  ymax = 3e4;
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
  ymax = 1e4;
end
[f] = plt_burst_diffets(ymax,thf03,thf04,thf05,hfinter03,hfinter04,hfinter05);
text(f.ax(1),0.06,0.5,strcat(num2str(round(perchf03)),'%'),'FontSize',12,'unit','normalized',...
    'horizontalalignment','left');
text(f.ax(2),0.06,0.5,strcat(num2str(round(perchf04)),'%'),'FontSize',12,'unit','normalized',...
    'horizontalalignment','left');
text(f.ax(3),0.06,0.5,strcat(num2str(round(perchf05)),'%'),'FontSize',12,'unit','normalized',...
    'horizontalalignment','left');
text(f.ax(1),0.04,0.9,'HF','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% print(f.fig,'-dpdf',strcat(figpath,'/selected_bursts_hf_',FLAG,'.pdf'));

%% get the refined catalog based on the tremor bursts
%%% depending on the what do you need, if you only need the time range, and in the next script, you
%%% just use the raw catalog, then add in the time range constraint and location constraint, we
%%% could stop here.
%%% But if you need also a refined catalog that is composed by the detections that fall into the
%%% time range of these bursts

if strcmp(TYPE, 'HF')  % opiton 1, use the HF bursts
  trange = [tranhf03;tranhf04;tranhf05];
end
  
hfnew = [];
lfnew = [];
for i = 1: size(trange, 1)
  indhf = find(hfbnd(:,daycol)==trange(i,1) & hfbnd(:,seccol)>=trange(i,2) & ...
               hfbnd(:,seccol)<=trange(i,3));
  hf = hfbnd(indhf,:);
  hfnew = [hfnew; hf];
  
  indlf = find(lfbnd(:,daycol)==trange(i,1) & lfbnd(:,seccol)>=trange(i,2) & ...
               lfbnd(:,seccol)<=trange(i,3));
  lf = lfbnd(indlf,:);
  lfnew = [lfnew; lf];
  
end
perchfnew = size(hfnew,1)/size(hfbnd,1)
perclfnew = size(lfnew,1)/size(lfbnd,1)


keyboard

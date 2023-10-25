% function [hfbnd, lfbnd] = catalog_of_region_of_interest(FLAG,TYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [hfbnd, lfbnd] = catalog_of_region_of_interest(FLAG,TYPE)
%
% -- This script is to read in the optional tremor catalog, either the raw
% merged catalog of all fams of interest, detected by PGC trio or LZB trio
% (FLAG) that contain false positives, or the refined catalog based on the
% temporally clustered tremor bursts from 'catalog_of_tremor_bursts.m'.
% -- Then it will further focus on a specific region of interest, and query
% the detections inside the study region. Those are 'repeating' detections
% that are active across several ETS episodes at the 'same' location.
% -- Analysing the spectra of these windows as a function of time would tell
% us, suppose the attenuation is the cause of the low corner-frequency (or
% low high-frequency energy), whether the attenuation is changing with time,
% ie., time-dependent.
% -- But this discussion is based on the assumption that the possible
% attenuation change in time could lead to a DETECTable change in corner-
% frequency of the spectra, which is supposed to be proportional to the
% reciprocal of the propagation time (propagation distance divided by the
% velocity).
% -- The region of interest could take the study of receiver function during
% ETS (Seismic evidence for megathrust fault-valve behavior during episodic
% tremor and slip, Gosselin_etal2020SciAdv).
% -- the resulting catalog not only is the catalog of detections, but also
% has a version of time windows defined by the earliest and lastest time of
% detections inside the bursts at the region of interest.
%
% Update 2021/10/26
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
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/28
% Last modified date:   2021/10/27
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
%     hffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxhf),'kmdcutnodou.burst',...
%       TYPE,'.',SUFFIXhf);
%     lffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxlf),'kmdcutnodou.burst',...
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
%     hffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',SUFFIXhf);
%     lffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',SUFFIXlf);
  hffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcut.',SUFFIXhf);
  lffname = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxlf),'kmdcut.',SUFFIXlf);

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

close all
   
%% define a region of interest
% All
famcol = 13;
hfplt = hftime(hftime(:,famcol)==2, :);
lfplt = lftime(lftime(:,famcol)==2, :);
xran = [-20 50];
yran = [-10 40];

[f4] = plt_cumulative_density_agu2021(hfplt,xran,yran);
print(f4.fig,'-dpdf',strcat(figpath,'/density_agu2021.pdf'));

[f4] = plt_cumulative_density(hfplt,lfplt,xran,yran);
% [f4] = plt_cumulative_density([],[],xran,yran);
% text(f4.ax(1), 0.85, 0.93, 'All','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);

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
% hold(f4.ax(1),'on');
% plot(f4.ax(1),bndse(:,1),bndse(:,2),'k-','linew',1.5);
% hold(f4.ax(1),'off');
% hold(f4.ax(2),'on');
% plot(f4.ax(2),bndse(:,1),bndse(:,2),'k-','linew',1.5);
% hold(f4.ax(2),'off');

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
% hold(f4.ax(1),'on');
% plot(f4.ax(1),bndcnt(:,1),bndcnt(:,2),'k-','linew',1.5);
% hold(f4.ax(1),'off');
% hold(f4.ax(2),'on');
% plot(f4.ax(2),bndcnt(:,1),bndcnt(:,2),'k-','linew',1.5);
% hold(f4.ax(2),'off');

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
% hold(f4.ax(1),'on');
% plot(f4.ax(1),bndnw(:,1),bndnw(:,2),'k-','linew',1.5);
% hold(f4.ax(1),'off');
% hold(f4.ax(2),'on');
% plot(f4.ax(2),bndnw(:,1),bndnw(:,2),'k-','linew',1.5);
% hold(f4.ax(2),'off');
       
% a rectangle N of the main box, close to fam 243
bndn = [-3 16;
       	-3 7;
        12 7;
        12 16;
        -3 16;
        ];
      
% hold(f4.ax(1),'on');
% plot(f4.ax(1),bndn(:,1),bndn(:,2),'r-','linew',2);
% hold(f4.ax(1),'off');
% hold(f4.ax(2),'on');
% plot(f4.ax(2),bndn(:,1),bndn(:,2),'r-','linew',2);
% hold(f4.ax(2),'off');

% a rectangle NE of the main box, close to fam 240
bndne = [12 18;
       	 12 8;
         18 8;
         18 18;
         12 18;
         ];
      
% hold(f4.ax(1),'on');
% plot(f4.ax(1),bndne(:,1),bndne(:,2),'c-','linew',2);
% hold(f4.ax(1),'off');
% hold(f4.ax(2),'on');
% plot(f4.ax(2),bndne(:,1),bndne(:,2),'c-','linew',2);
% hold(f4.ax(2),'off');

% plot station names
dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');

stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4];
lzbtrio = stainfo(ind,:);
for i = 1: 3
    [dxlzb(i),dylzb(i)] = absloc2relaloc(str2num(lzbtrio(i,3)),str2num(lzbtrio(i,2)),...
                    -123.585000, 48.436667);
end

ind=[5,8,6];
pgctrio = stainfo(ind,:);
for i = 1: 3
    [dxpgc(i),dypgc(i)] = absloc2relaloc(str2num(pgctrio(i,3)),str2num(pgctrio(i,2)),...
                    -123.585000, 48.436667);
end

otherind = setdiff(1:12, union([3,11,4],[5,8,6]));
othersta = stainfo(otherind,:);
for i = 1: length(otherind)
    [dxother(i),dyother(i)] = absloc2relaloc(str2num(othersta(i,3)),str2num(othersta(i,2)),...
                    -123.585000, 48.436667);
end

%plot stations and names
ax = f4.ax(1);
hold(ax,'on');
for i=1:3
    % p1 is lzb trio stations
    p1 = plot(ax,dxlzb(i),dylzb(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','k');
    text(ax,dxlzb(i)+1,dylzb(i)+1,lzbtrio(i,1));
    % p2 is pgc trio stations
    p2 = plot(ax,dxpgc(i),dypgc(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','r');
    text(ax,dxpgc(i)+1,dypgc(i)+1,pgctrio(i,1)); 
end

for i = 1: length(otherind)
  % p3 is other unused stations
    p3 = plot(ax,dxother(i),dyother(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','w','clipping','on');
    text(ax,dxother(i)+1,dyother(i)+1,othersta(i,1),'clipping','on');
end

%plot fam names
loccont = [
           -123.585000 48.436667 36.8800;   % 002   
           -123.549500 48.540833 38.5600;   % 243
           -123.382333 48.574167 40.9800;   % 240
           ];
[relacont(:,1), relacont(:,2)] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.585000, 48.436667);
scatter(ax,relacont(:,1), relacont(:,2),30,'k','filled');
for i = 1: size(loccont,1)
  text(ax,relacont(i,1)+1, relacont(i,2)+1, nfampool(i,:),'fontsize',8,'backgroundColor','w','Margin',0.5);
end

hold(ax,'off');

print(f4.fig,'-dpdf',strcat(figpath,'/density',FLAG,'_',TYPE,'_ETSallraw.pdf'));

% keyboard


%% output the catalog of detections inside the region and time ranges
bnduse = bndnw;
if isequal(bnduse,bnd)
  bndflag = 'bnd';
elseif isequal(bnduse,bndse)
  bndflag = 'bndse';
elseif isequal(bnduse,bndcnt)
  bndflag = 'bndcnt';
elseif isequal(bnduse,bndnw)
  bndflag = 'bndnw';
elseif isequal(bnduse,bndn)
  bndflag = 'bndn';
elseif isequal(bnduse,bndne)
  bndflag = 'bndne';
end

% step 1, events further inside the bound
[iin,ion] = inpolygon(hftime(:,1),hftime(:,2),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
hfbnd = hftime(isinbnd == 1, :);

[iin,ion] = inpolygon(lftime(:,1),lftime(:,2),bnduse(:,1),bnduse(:,2));
isinbnd = iin | ion;
lfbnd = lftime(isinbnd == 1, :);

% All
hfplt = hfbnd;
lfplt = lfbnd;
xran = [-25 30];
yran = [-20 30];
[f4] = plt_cumulative_density(hfplt,lfplt,xran,yran);
text(f4.ax(1), 0.85, 0.93, strcat({'All in '},bndflag),'FontSize',12,'unit','normalized',...
  'horizontalalignment','center','EdgeColor','k','Margin',2);

% %%%%%%%%%%%% OLD parts when running 'catalog_of_tremor_bursts' %%%%%%%%%%%%%%%%%
% hfnew = [];
% lfnew = [];
% ntrange = [];
% k = 0;
% for i = 1: size(trange,1)  
%   % step 1, events inside the range of bursts from the entire catalog
%   indhf = find(hftime(:,daycol)==trange(i,1) & hftime(:,seccol)>=trange(i,2) & ...
%                hftime(:,seccol)<=trange(i,3));
%   hf = hftime(indhf,:);
%   
%   indlf = find(lftime(:,daycol)==trange(i,1) & lftime(:,seccol)>=trange(i,2) & ...
%                lftime(:,seccol)<=trange(i,3));
%   lf = lftime(indlf,:);
%   
%   % step 2, events further inside the bound
%   [iin,ion] = inpolygon(hf(:,1),hf(:,2),bnduse(:,1),bnduse(:,2));
%   isinbnd = iin | ion;
%   hfbnd = hf(isinbnd == 1, :);
%   
%   [iin,ion] = inpolygon(lf(:,1),lf(:,2),bnduse(:,1),bnduse(:,2));
%   isinbnd = iin | ion;
%   lfbnd = lf(isinbnd == 1, :);
%     
%   % new catalogs
%   hfnew = [hfnew; hfbnd];
%   lfnew = [lfnew; lfbnd];
%   
% %   % new time ranges, should be shorter than the old ones individually,
% %   % use HF as reference only
% %   if ~isempty(hfbnd)
% %     k = k+1;
% %     ntrange(k,1) = trange(i,1);
% %     ntrange(k,2) = min(hfbnd(:,seccol));
% %     ntrange(k,3) = max(hfbnd(:,seccol));
% %   end
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%% maybe a newer scheme for the start and end time of the new time range?
%   %%% there will be some extreme outliers dominating the range by using max
%   %%% try 'prctile'
%   if ~isempty(hfbnd)
%     k = k+1;
%     ntrange(k,1) = trange(i,1);
%     ntrange(k,2) = prctile(hfbnd(:,seccol),2);
%     ntrange(k,3) = prctile(hfbnd(:,seccol),98);
%   end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% end
% %%%%%%%%%%%% OLD parts when running 'catalog_of_tremor_bursts' %%%%%%%%%%%%%%%%%


if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  hffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxhf),'dcutnodb.',...
    TYPE,bndflag,'.',SUFFIXhf);
  lffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxlf),'dcutnodb.',...
    TYPE,bndflag,'.',SUFFIXlf);  
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
  hffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxhf),'dcutnodb.',...
    TYPE,bndflag,'.',SUFFIXhf);
  lffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxlf),'dcutnodb.',...
    TYPE,bndflag,'.',SUFFIXlf);
end

fidhf = fopen(hffname,'w+');
fprintf(fidhf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        hfbnd');
fclose(fidhf);

fidlf = fopen(lffname,'w+');
fprintf(fidlf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        lfbnd');
fclose(fidlf);

% keyboard

  







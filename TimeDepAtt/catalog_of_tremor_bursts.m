function [trange, hfnew, lfnew] = catalog_of_tremor_bursts(FLAG,TYPE,bndflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [trange, hfnew, lfnew] = catalog_of_tremor_bursts(FLAG,TYPE,bndflag)
% -- This script is to read in the merged catalog of all fams of interest, 
% detected by PGC trio or LZB trio (FLAG), then use the same time-clustering
% algorithm in the 'identify_tremor_bursts_intertime.m' to obtain all the 
% temporally clustered tremor bursts. Note that it does NOT guarantee the 
% spatial clustering of the detections. 
% -- The ranges of bursts are obtained separately from HF and LF (optional), in
% case that LF detected some bursts with better SNR that are not detected in 
% HF.
% -- The detections inside the time ranges defined by these bursts then become
% the refined catalog, which are going to be the other option of the start of
% our time-dependent attenuation study. From this refined catalog, we will 
% focus on the entire region imaged, or a small region, to analyze the spectra
% of tremor windows whose location fall into that region of interest in each ETS.
% -- Note that if the region of interest is smaller than imaged, the catalog
% has to be truncated, meaning that the times are narrower. After all, I think
% we care less about the detections, but the time ranges that will be divided
% into consecutive windows to analyse the spectra.
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
%
% Note 2021/10/27:
% --Using MAD for the time separation threshold for the same data set seems fine,
%   but it is problematic when comparing diff data sets. A sparser set would
%   lead to a higher threshold, using which would result in more and longer
%   burst windows that is unfair for a denser and more clustered set. Example
%   could be made between 'bnd' and 'bndn' where 'bnd' clearly has more and more
%   clustered detections. Using the same 8*MAD of themselves lead to a unfairly
%   longer windows for 'bndn', this makes no sense to me. It looks like using
%   the same constant threshold for all subregions is better.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/26
% Last modified date:   2021/10/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
% clear
close all

%% default value for easy debugging
defval('FLAG','PGC');
defval('TYPE','HF');
defval('bndflag','bnd');

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
    hffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxhf),'dcutnodb.',...
      TYPE,bndflag,'.',SUFFIXhf);
    lffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxlf),'dcutnodb.',...
      TYPE,bndflag,'.',SUFFIXlf);  

    
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
    hffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxhf),'dcutnodb.',...
      TYPE,bndflag,'.',SUFFIXhf);
    lffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxlf),'dcutnodb.',...
      TYPE,bndflag,'.',SUFFIXlf);

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

% obtain the relative time in days
ih03 = find(hftime(:,daycol) < 2004*1000);
hftime(ih03,ncol+1) = (hftime(ih03,daycol)-2003060)+hftime(ih03,seccol)./(3600.*24);
ih04 = find(hftime(:,daycol) < 2005*1000 & hftime(:,daycol) > 2004*1000);
hftime(ih04,ncol+1) = (hftime(ih04,daycol)-2004194)+hftime(ih04,seccol)./(3600.*24);
ih05 = find(hftime(:,daycol) > 2005*1000);
hftime(ih05,ncol+1) = (hftime(ih05,daycol)-2005254)+hftime(ih05,seccol)./(3600.*24);
% % rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% % to down-dip and strike
% [hftime(ih05,1),hftime(ih05,2)] = coordinate_rot(hftime(ih05,1),hftime(ih05,2),rotang,0,0);

il03 = find(lftime(:,daycol) < 2004*1000);
lftime(il03,ncol+1) = (lftime(il03,daycol)-2003060)+lftime(il03,seccol)./(3600.*24);
il04 = find(lftime(:,daycol) < 2005*1000 & lftime(:,daycol) > 2004*1000);
lftime(il04,ncol+1) = (lftime(il04,daycol)-2004194)+lftime(il04,seccol)./(3600.*24);
il05 = find(lftime(:,daycol) > 2005*1000);
lftime(il05,ncol+1) = (lftime(il05,daycol)-2005254)+lftime(il05,seccol)./(3600.*24);
% % rotate to the strike N30W and dip N60E direction,now the 1,2 col is changing from E(043) and N(043),
% % to down-dip and strike
% [lftime(il05,1),lftime(il05,2)] = coordinate_rot(lftime(il05,1),lftime(il05,2),rotang,0,0);

% keyboard

%% separation in time between itself and its preceding
% for the detections, obtain the separation in time between itself and its preceding
% detection for the catalog of each ETS.
% 1st col: occurence time
% 2nd col: inter-event time to its preceding detection
hfinter03 = interevt_time(hftime(ih03,end));
hfinter04 = interevt_time(hftime(ih04,end));
hfinter05 = interevt_time(hftime(ih05,end));
hfinter = [hfinter03; hfinter04; hfinter05];

lfinter03 = interevt_time(lftime(il03,end));
lfinter04 = interevt_time(lftime(il04,end));
lfinter05 = interevt_time(lftime(il05,end));
lfinter = [lfinter03; lfinter04; lfinter05];

%% plot the inter-detection time
%%% HF catalog
% set a threshold of inter-detection time
ttol1 = 1e-3*ones(3,1); 
ttol2 = 2e-3*ones(3,1);
% ttol2 = 8*[mad(hfinter03(:,2), 1); mad(hfinter04(:,2), 1); mad(hfinter05(:,2), 1)];
ymax = 1.2*max([ttol2; ttol1]);
[f] = plt_interevt_time_diffets(hfinter03,hfinter04,hfinter05,ymax,ttol1,ttol2);
text(f.ax(1),0.04,0.9,'HF','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
print(f.fig,'-dpdf',strcat(figpath,'/inter-detection_time_HF',FLAG,bndflag,'.pdf'));  


% %%% LF catalog
% ttol1 = 3e-3*ones(3,1);
% % ttol2 = 2e-3;
% ttol2 = 8*[mad(lfinter03(:,2), 1); mad(lfinter04(:,2), 1); mad(lfinter05(:,2), 1)];
% ymax = 1.2*max([ttol2; ttol1]);
% [f] = plt_interevt_time(ttol1,ttol2,ymax,lfinter03,lfinter04,lfinter05);
% text(f.ax(1),0.04,0.9,'LF','FontSize',11,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% print(f.fig,'-dpdf',strcat(figpath,'/inter-detection_time_LF',FLAG,bndflag,'.pdf')); 

% keyboard

%% plot the detections in map to see if they are clustered in space
%%% this is to make sure the threshold of separation time indeed capture most of the clustered
%%% bursts, as this time-clustering algorithm does not guarantee a cluster in space at the same time
CHECK = 0;
% CHECK = 1;

if CHECK
  checkran = [
%               2003061 86400*0.25 86400*0.35;
              2003061 86400*0.28 86400*0.35;
              2003061 39714.0000000000 40593.0000000000;
              2003061 86400*0.48 86400*0.52;
              2003061 86400*0.53 86400*0.57;
              2003061 86400*0.58 86400*0.63;
              2003061 52043.0000000000 52572.0000000000;
              ];
  
  for i = 1: size(checkran,1)
    disp(checkran(i,:));
    colnum = [daycol seccol 1 2]; 
    
    
    %%% define and position the figure frame and axes of each plot
    f.fig=figure;
    widin = 5;  % maximum width allowed is 8.5 inches
    htin = 5;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    nrow = 1;
    ncol = 1;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end
    
    [f.ax(1)] = plt_detections_inmap(f.ax(1),hftime,checkran(i,:),colnum);
  end
  
  
end

%% group the tremor bursts 
%%% group the indice of the tremor bursts accroding to the threshold on the inter-detection time and
%%% number of detections in the burst for both HF and LF
%%%% HF 
if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  % set a threshold of inter-detection time
  ttolhf = 1e-3;
%   ttolhf = 5e-4;
  % set a threshold of minimum number of detections in the burst period
%   ntolhf = 15;
%   ntollf = 15;
  ntolhf = 30;
  ntollf = 15;
  [bursthf03, nburstlf03] = group_tremor_burst(hfinter03,lftime(il03,end),ttolhf,ntolhf,0);
  [bursthf04, nburstlf04] = group_tremor_burst(hfinter04,lftime(il04,end),ttolhf,ntolhf,0);
  [bursthf05, nburstlf05] = group_tremor_burst(hfinter05,lftime(il05,end),ttolhf,ntolhf,0);
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
%   ttolhf = 1e-3;  % param used until 2021/10/24
%   ntolhf = 30;  % param used until 2021/10/24
  ttolhf = 2e-3*ones(3,1);  % param used from 2021/10/27
%   ttolhf = 6*[mad(hfinter03(:,2), 1); mad(hfinter04(:,2), 1); mad(hfinter05(:,2), 1)]; 
  if strcmp(bndflag,'bndse') || strcmp(bndflag,'bndcnt') || strcmp(bndflag,'bndnw')
    ntolhf = 8;  % param used from 2021/10/27, these are much smaller regions
  else
    ntolhf = 15;  % param used from 2021/10/27
  end
  ntollf = 10;
  [bursthf03, nburstlf03] = group_tremor_burst(hfinter03,lftime(il03,end),ttolhf(1),ntolhf,0);
  [bursthf04, nburstlf04] = group_tremor_burst(hfinter04,lftime(il04,end),ttolhf(2),ntolhf,0);
  [bursthf05, nburstlf05] = group_tremor_burst(hfinter05,lftime(il05,end),ttolhf(3),ntolhf,0);
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
perclf = (nlf03+nlf04+nlf05)/(length(il03)+length(il04)+length(il05))*100;

% keyboard

% %%%% LF 
% if strcmp(FLAG, 'TWKB')     % use TWKB catalog
%   ttollf = 2e-3;
% elseif strcmp(FLAG, 'PGC')     % use PGC catalog
%   ttollf = 2e-3;
% end
% [burstlf03, nbursthf03] = group_tremor_burst(lfinter03,hftime(ih03,end),ttollf,ntollf,0);
% [burstlf04, nbursthf04] = group_tremor_burst(lfinter04,hftime(ih04,end),ttollf,ntollf,0);
% [burstlf05, nbursthf05] = group_tremor_burst(lfinter05,hftime(ih05,end),ttollf,ntollf,0);
% if ~isempty(burstlf03)
%   [tlf03,tranlf03,nlf03] = burst_range(burstlf03,lfinter03,2003060);
% else
%   tlf03 = []; tranlf03 = []; nlf03 = 0;
% end
% if ~isempty(burstlf04)
%   [tlf04,tranlf04,nlf04] = burst_range(burstlf04,lfinter04,2004194);
% else
%   tlf04 = []; tranlf04 = []; nlf04 = 0;
% end
% if ~isempty(burstlf05)
%   [tlf05,tranlf05,nlf05] = burst_range(burstlf05,lfinter05,2005254);
% else
%   tlf05 = []; tranlf05 = []; nlf05 = 0;
% end

% keyboard


%% plot the selected tremor bursts
%%%% HF
if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  ymax = 3e4;
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
  ymax = 4e3;
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
print(f.fig,'-dpdf',strcat(figpath,'/selected_bursts_HF',FLAG,bndflag,'.pdf'));

% %%%% LF
% if strcmp(FLAG, 'TWKB')     % use TWKB catalog
%   ymax = 3e4;
% elseif strcmp(FLAG, 'PGC')     % use PGC catalog
%   ymax = 1e3;
% end
% [f] = plt_burst_diffets(ymax,tlf03,tlf04,tlf05,lfinter03,lfinter04,lfinter05);
% % text(f.ax(1),0.06,0.8,strcat(num2str(round(perclf03)),'%'),'FontSize',12,'unit','normalized',...
% %     'horizontalalignment','left');
% % text(f.ax(2),0.06,0.8,strcat(num2str(round(perclf04)),'%'),'FontSize',12,'unit','normalized',...
% %     'horizontalalignment','left');
% % text(f.ax(3),0.06,0.8,strcat(num2str(round(perclf05)),'%'),'FontSize',12,'unit','normalized',...
% %     'horizontalalignment','left');
% text(f.ax(1),0.04,0.9,'LF','FontSize',11,'unit','normalized','horizontalalignment','left',...
%     'EdgeColor','k','Margin',2);
% print(f.fig,'-dpdf',strcat(figpath,'/selected_bursts_LF',FLAG,bndflag,'.pdf'));

% keyboard

% %% get a union set of the bursts imaged by either HF and LF
% %%% allowing some difference in time of start and end of the same burst, as they are unlikely to be
% %%% exactly the same
% %%% NOTE that there is a flaw here, and it may not be the problem of 'merge_burst' function.
% %%% Although the bursts grouped separately by hf and lf are themselves nonpverlapping in time, but
% %%% when you combine them, the union set then could overlap, and this function would merge the
% %%% overlapping bursts from either hf or lf into a longer burst
% %%% So this would leave us 2 options,
% %%% option 1:  regard the hf bursts as reliable, and use them only
% %%% option 2:  use the merged bursts which are individually longer in time, but smaller in total number
% [ntran03,merind03,mertran03] = merge_burst(tranhf03,tranlf03);
% [ntran04,merind04,mertran04] = merge_burst(tranhf04,tranlf04);
% [ntran05,merind05,mertran05] = merge_burst(tranhf05,tranlf05);
% 
% % [f] = plt_burst(tlf03,tlf04,tlf05,lfinter03,lfinter04,lfinter05);
% % text(f.ax(1),0.04,0.9,'HF','FontSize',11,'unit','normalized','horizontalalignment','left',...
% %     'EdgeColor','k','Margin',2);
% % print(f.fig,'-dpdf',strcat(figpath,'/merged_bursts_hf_',FLAG,'.pdf'));
% 
% % keyboard

%% get the refined catalog based on the tremor bursts
%%% depending on the what do you need, if you only need the time range, and in the next script, you
%%% just use the raw catalog, then add in the time range constraint and location constraint, we
%%% could stop here.
%%% But if you need also a refined catalog that is composed by the detections that fall into the
%%% time range of these bursts

if strcmp(TYPE, 'HF')  % opiton 1, use the HF bursts
  trange = [tranhf03;tranhf04;tranhf05];
elseif strcmp(TYPE, 'LF')   % opiton 2, use the LF bursts
  trange = [tranlf03;tranlf04;tranlf05];
elseif strcmp(TYPE, 'Common')   % opiton 3, use the merged bursts common to both HF and LF
  trange = [ntran03;ntran04;ntran05];
end

% get the new catalog composed only by the bursts
hfnew = refine_catalog_by_trange(hftime,daycol,seccol,trange);
lfnew = refine_catalog_by_trange(lftime,daycol,seccol,trange);

perchfnew = size(hfnew,1)/size(hftime,1)
perclfnew = size(lfnew,1)/size(lftime,1)

if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  hffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxhf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXhf);
  lffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxlf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXlf);  
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
  hffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxhf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXhf);
  lffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxlf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXlf);
end

fidhf = fopen(hffname,'w+');
fprintf(fidhf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f %.5f \n',...
        hfnew');
fclose(fidhf);

fidlf = fopen(lffname,'w+');
fprintf(fidlf,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f %.5f \n',...
        lfnew');
fclose(fidlf);

fid = fopen(strcat(hypopath, '/tremor_burst_ranges_',FLAG,TYPE,bndflag),'w+');
fprintf(fid,'%d %9.2f %9.2f \n',trange');
fclose(fid);

keyboard











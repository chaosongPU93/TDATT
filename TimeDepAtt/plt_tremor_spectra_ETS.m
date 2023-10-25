% function plt_tremor_spectra_ETS(FLAG,TYPE,fam,bndflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plt_tremor_spectra_ETS(FLAG,TYPE,fam,bndflag)
% -- This script is to read in the tremor catalog descripted as the windows
% of time, resulting from 'catalog_of_region_of_interest.m', so that we know
% which portion of daily tremor data to analyse.
% -- The windows denote the times when the region of interest is active in 
% different ETS episodes. Then we just sort the spectra in time.
% -- All the windows in one ETS, however you bin them, represent one 
% realization of our experiment. The binning method can be optional, either
% by each tremor burst, or by an equal length of time (??). But note that
% each burst has a variable length, so its spectra is an 'average' of several
% equal-length subwindows. Then we just plot each of the bursts in the same
% subplot and extract an 'averge' spectra from all of them. 
% -- but I strongly suspect there will be the same bursts occurring every ETS
% so that you could find an exactly one-to-one correspondence. So maybe obtain
% an average is a better idea.
% -- Q: however, what if the same location has amplitude high or low varying
% with time? Then get an direct average or median is less meaningful. I guess
% we'd better normalize the amplitude within a frequency range of interest, 
% then we could wither look at the each spectrum plotting on top of each other,
% or compute an average or median from all normalized spectra.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/28
% Last modified date:   2021/10/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
clear
close all

%% default value for easy debugging
defval('FLAG','PGC');
defval('TYPE','HF');
defval('fam','002');
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
%   if strcmp(CATA, 'raw')	% option 1: the raw merged catalog
%     hffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',SUFFIXhf);
%     lffname = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',SUFFIXlf);
%   elseif strcmp(CATA, 'burst')  % option 2: the refined catalog based on tremor bursts
  hffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxhf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXhf);
  lffname = strcat(hypopath, '/eloc.lzb.pj21.',num2str(distmaxlf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXlf);  
%   end
  
  % read the rotation parameters accroding to the fam
  % Note: here I guess we have to choose a controlling fam and use its rotation para to read the data,
  % etc.
  % get permanent and polaris station rotation parameters
  sft2=0;     % centroid shift of station 2
  sft3=0;     % centroid shift of station 3
  CATA = 'new';
  [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
  reftime = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;
  
  PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
  POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

  
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
    'SILB '
    'KLNB '
    ];
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
%   elseif strcmp(CATA, 'burst')  % option 2: the refined catalog based on tremor bursts
  hffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxhf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXhf);
  lffname = strcat(hypopath, '/eloc.pgc.pj21.',num2str(distmaxlf),'dcutnodb.bst',...
    TYPE,bndflag,'.',SUFFIXlf);
%   end
  
  % read the rotation parameters accroding to the fam
  % Note: here I guess we have to choose a controlling fam and use its rotation para to read the data,
  % etc.
  % get permanent and polaris station rotation parameters
  if isequal(fam,'002')
    %     [timoffrot,stas,bostname,tempoffs] = GetDays(fam,freqflag);
    %     nsta=size(stas,1);         %  number of stations
    % get permanent and polaris station rotation parameters
    sft2=0;     % centroid shift of station 2
    sft3=0;     % centroid shift of station 3
    CATA = 'fixed';
    [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
  else
    sft2=0;     % centroid shift of station 2
    sft3=0;     % centroid shift of station 3
    CATA = 'new';
    [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
    reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
    PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
    POLROTS(:,4) = POLROTS(:,4)-reftime;
  end  
  PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;     % convert 2-3 columns to rad
  POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;

end
% hftime = load(hffname);
% lftime = load(lffname);

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
  
% % sort according to day, sec of the window center
% daycol = 14;
% seccol = 15;
% hftime = sortrows(hftime, [daycol, seccol]);
% lftime = sortrows(lftime, [daycol, seccol]);
% 
% ncol = size(hftime,2);

%%% load the time range of the tremor bursts
%%% option 1, windows from automatic clustering algorithm
fname = strcat(hypopath, '/tremor_burst_ranges_',FLAG,TYPE,bndflag);
trange = load(fname);

% %%% option 2, windows from Allan's eyes using ginput
% jdays= [2003 062; %jdays are for data.
%         2003 063;
%         2004 196;
%         2004 197;
%         2005 254;
%         ];
% ntrange = [];      
% for id = 1: size(jdays,1)        
%   wins=getspecwins(jdays(id, :));
%   ran = zeros(size(wins,1), 3);
%   ran(:, 1) = (jdays(id, 1)*1000 + jdays(id, 2))*ones(size(wins,1),1);
%   ran(:, 2:3) = wins;
%   ntrange = [ntrange; ran];
% end


%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];

% keyboard

%% For each burst in ETS, plot the spectra
% for filtering
sps = 100;
hi=19;
lo=.1; 
npo=2;
npa=2;

% dates and ets
%%% NOTE: 'dates' is the dates that the tremor is active at the region of interest, so that the
%%% waveform at the stations on these dates you are looking at mainly comes from the sources at the 
%%% region of interest. We want to see if there is a noticable change in spectra during the burst 
%%% windows on these dates
dates = unique(trange(:,1));
years = unique(floor(dates/1000));
nets = length(years);

% keyboard
for iets = 1: nets
  % dates in each ets
  year = years(iets);
  datesets = dates(floor(dates/1000)==year);
  
  if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
    stas(4,:)='KELB ';
  else
    stas(4,:)='KLNB ';  % remember to change it back
  end
  
  k = 0;  % for counting
  optcatdat = [];   % concataneted array for all windows in the same ETS, optimal component
  ortcatdat = [];   % orthogonal component
  optsnrcatdat = [];  % optimal component within a freq range of a high SNR 
  for i = 1: length(datesets)
    
    %%% read daily tremor data
    date = datesets(i);
    jday = floor(date-year*1000);
    
    [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
                              PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);
%     %data within a good SNR freq range 
%     [STAoptsnr,~,~,~] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
%                               PERMROTS,POLROTS,sps,2.5,5,npo,npa,[],[],[],[]);
    if fileflag == 0    % means there are missing files
      fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
      continue    % continue to the next day
    end
                                
    %%% We do 2 things here:
    %%% 1. simply get the spectra of all windows and store them
    %%% 2. concatenate all windows from the same ETS as if they are continuous, then get the
    %%% spectrogram of it
    % get the windows (time ranges) of the same day, for efficiency
    rangetemp = trange(trange(:,1)==datesets(i), :);
%     keyboard
    for j = 1: size(rangetemp,1)
      tst = rangetemp(j,2);
      ted = rangetemp(j,3);
      wlen = 2048;
      wlensec = wlen/sps;
      if ted-tst>=wlensec
        polap = 50;
        nfft = wlen;
        Fs = sps;
        optdat = STAopt(max(floor(tst*sps),1): min(floor(ted*sps-1),86400*sps), :);
        ortdat = STAort(max(floor(tst*sps),1): min(floor(ted*sps-1),86400*sps), :);
%         optsnrdat = STAoptsnr(max(floor(tst*sps),1): min(floor(ted*sps-1),86400*sps), :);
        for ista = 1: nsta
          [pcopt(:,ista),pcft]=pchave(optdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
%           [pcort(:,ista),~]=pchave(ortdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
%           [pcoptsnr(:,ista),~]=pchave(optsnrdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
        end
        
        % stores for every rtm
        k = k+1;
        pcetsopt(:,:,k) = pcopt(:,:);   % ifreq, ista, irtm
%         pcetsort(:,:,k) = pcort(:,:);
%         pcetsoptsnr(:,:,k) = pcoptsnr(:,:);
        
        % concatenate all windows in the same ETS
        optcatdat = [optcatdat; optdat];
%         ortcatdat = [ortcatdat; ortdat];
%         optsnrcatdat = [optsnrcatdat; optsnrdat];

      end  % if length is enough
%       keyboard
    end % bursts in each date
    
  end % dates in each ets
  
  % store spectra of all windows of the different ETSs separately
  pcall{iets,1} = pcetsopt(:,:,:);
%   pcall{iets,2} = pcetsort(:,:,:);
%   pcall{iets,3} = pcetsoptsnr(:,:,:);
  
  % now for the concatenated window, get the spectrogram as if it is continuous
  seglen = wlen;
  %         nfft = pow2(nextpow2(seglen));
  nfft = wlen;
  window = hann(nfft);
  Fs = sps;
  olap = polap/100;
  nolap = seglen*olap;
  
  spgcatopt = [];
  spgcatort = [];
  spgcatoptsnr = [];
  for ista = 1: nsta
    [spgcatopt(:,:,ista), spgft, t, psddB]=spectrogram2(optcatdat(:,ista+1),nfft,Fs,seglen,nolap,'s');
%     [spgcatort(:,:,ista), ~,~,~]=spectrogram2(ortcatdat(:,ista+1),nfft,Fs,seglen,nolap,'s');
%     [spgcatoptsnr(:,:,ista), ~,~,~]=spectrogram2(optsnrcatdat(:,ista+1),nfft,Fs,seglen,nolap,'s');
  end
  
  % store the spectrograms of the different ETSs separately
  spgall{iets,1} = spgcatopt(:,:,:);
%   spgall{iets,2} = spgcatort(:,:,:);
%   spgall{iets,3} = spgcatoptsnr(:,:,:);
  tall{iets,1} = t;
  
end % all ets

% clean some temp variables to make more room
clear STAopt STAort 
clear optdat ortdat 
clear pcopt pcort 
clear pcetsopt pcetsort 
clear spgcatopt spgcatort 

keyboard


%% plot the spectra of optimal components of all ETS
xran = [0.1 20];
yran = [1e-2 1e2];
minfnorm = 1.5; maxfnorm = 5;
[f] = plt_spectra_of_bursts_norm(years,stas,pcft,pcall,minfnorm,maxfnorm,xran,yran);
%picked corner frequency from the median spectra
fc_obs = [4.54; 4.0; 4.0; 5.27];

% xvect = [fc_obs(1) fc_obs(1)];
% yvect = [2 10];
% drawArrow(f.ax(9),xvect,yvect,xran,yran,'linewidth',1);
plot(f.ax(9),[fc_obs(1) fc_obs(1)],yran,'c--','linew',1.5);
plot(f.ax(10),[fc_obs(2) fc_obs(2)],yran,'c--','linew',1.5);
plot(f.ax(11),[fc_obs(3) fc_obs(3)],yran,'c--','linew',1.5);
plot(f.ax(12),[fc_obs(4) fc_obs(4)],yran,'c--','linew',1.5);

print(f.fig,'-dpdf',strcat(figpath,'/aaaa_spa',FLAG,bndflag,'.pdf'));
% print(f.fig,'-dpdf',strcat(figpath,'/bbbb_spectra',FLAG,bndflag,'.pdf'));

% keyboard


%% plot the spectrogram of concatanated windows for all ETS 
for iets = 1: nets

  yran = [2 6];
  minfnorm = 1.5; maxfnorm = 5;
  filtflag = 1;
  if filtflag
    cran = [-0.2 0.2];
    filtsigma = [2 10];
  else
    cran = [-0.4 0.4];
    filtsigma = [];
  end
  
  if iets==1   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
    stas(4,:)='KELB ';
  else
    stas(4,:)='KLNB ';  % remember to change it back
  end
  [f] = plt_spectrogram_of_bursts_norm(years,stas,iets,spgft,spgall,tall,minfnorm,maxfnorm,yran,cran,...
                                              filtflag,filtsigma);  
  plot(f.ax(1),f.ax(1).XLim, [4.8 4.8],'k--');
  plot(f.ax(2),f.ax(2).XLim, [2.55 2.55],'k--');
  plot(f.ax(2),f.ax(2).XLim, [4.4 4.4],'k--');
  plot(f.ax(3),f.ax(3).XLim, [2.55 2.55],'k--');
  plot(f.ax(3),f.ax(3).XLim, [3.4 3.4],'k--');
  plot(f.ax(3),f.ax(3).XLim, [4.4 4.4],'k--');
  plot(f.ax(4),f.ax(4).XLim, [2.9 2.9],'k--');
  plot(f.ax(4),f.ax(4).XLim, [4.0 4.0],'k--');
  
  if filtflag
    str = 'spg_gauf';
  else
    str = 'spg';
  end
  print(f.fig,'-dpdf',strcat(figpath,'/aaaa_',str,num2str(iets),FLAG,bndflag,'.pdf'));
%   print(f.fig,'-dpdf',strcat(figpath,'/bbbb_spectrogram',num2str(iets),FLAG,bndflag,'.pdf'));
  
end   % all ets

% keyboard
% close all

%% focus on the band feature at SILB, zoom-in; strech vertically; reference line; ticks on right
%%% it seems that the band feature is very obvious at stations SILB (2 bands), KLNB (1 band), and
%%% SSIB (1 band). At SILB, especially there seems to be a positive slope in the band, indicating a
%%% drifting in corner frequency (increasing), which may correspond to a increase in the velocity in
%%% the low-velocity layer in the upper-most crust of the subducting plate that is right beneath the
%%% tremor region (interface). The positive slope seems to exist at KLNB too, but less apparent at
%%% SSIB.
for iets = 1: nets
  for ista = 3: 3
    %   ista = 4;
    
    yran = [2 6];
    minfnorm = 1.5; maxfnorm = 5;
    filtflag = 1;
    if filtflag
      cran = [-0.15 0.15];
      % control the smoothing degree in each dimension, larger sigma mean higher smoothing; the size
      % of the filter as default is 2*ceil(2.*sigma)+1. Using the same sigma but increase size is 
      % useless as the gaussian distribution is nearly zero beyond 3*sigma 
      filtsigma = [2 15];
    else
      cran = [-0.3 0.3];
      filtsigma = [];
    end

    [f,tets,fgmax,fgmaxfit,fran] = plt_spgram_of_bursts_norm_zoom(years,stas,iets,ista,spgft,spgall,...
                                         tall,minfnorm,maxfnorm,yran,cran,filtflag,filtsigma);
    if ~filtflag
      ax = f.ax(1);
      hold(ax,'on');
      if ista == 1
        plot3(ax,ax.XLim, [4.8 4.8],10*ones(2,1),'k--');
      elseif ista == 2
        plot3(ax,ax.XLim, [2.55 2.55],10*ones(2,1),'k--');
        plot3(ax,ax.XLim, [4.4 4.4],10*ones(2,1),'k--');
%         plot(ax,ax.XLim, [3 3],10*ones(2,1),'w--');
%         plot(ax,ax.XLim, [3 3.1],10*ones(2,1),'w--');
      elseif ista == 3
        plot3(ax,ax.XLim, [2.55 2.55],10*ones(2,1),'k--');
        plot3(ax,ax.XLim, [3.4 3.4],10*ones(2,1),'k--');
        plot3(ax,ax.XLim, [4.4 4.4],10*ones(2,1),'k--');
%         plot3(ax,ax.XLim, [2.8 2.8],10*ones(2,1),'w--');
%         plot3(ax,ax.XLim, [2.8 3],10*ones(2,1),'w--');
%         plot3(ax,ax.XLim, [3.7 3.7],10*ones(2,1),'w--');
%         plot3(ax,ax.XLim, [3.7 4.0],10*ones(2,1),'w--');
      elseif ista == 4
        plot3(ax,ax.XLim, [2.9 2.9],10*ones(2,1),'k--');
        plot3(ax,ax.XLim, [4.0 4.0],10*ones(2,1),'k--');
%         plot3(ax,ax.XLim, [3.5 3.5],10*ones(2,1),'w--');
%         plot3(ax,ax.XLim, [3.5 3.7],10*ones(2,1),'w--');
      end
    end
    
    if filtflag
      str = 'spg_gauf';
    else
      str = 'spg';
    end
    print(f.fig,'-dpdf',strcat(figpath,'/aaaa_',str,'_zoom',num2str(iets),...
      strtrim(stas(ista, :)),FLAG,bndflag,'.pdf'));
%   print(f.fig,'-dpdf',strcat(figpath,'/bbbb_spectrogram',num2str(iets),FLAG,bndflag,'.pdf'));
%   print(f.fig,'-dpdf',strcat(figpath,'/aa',num2str(iets),'.pdf'));

    if filtflag    
      yran = [3 4.5];
      [f,fitobj,~,~,fgmaxfit] = plt_spgram_fcfit(tets,fgmax(:,1),fran(:,1:2),yran);
      text(0.04,0.95,num2str(years(iets)),'FontSize',10,'unit','normalized',...
        'horizontalalignment','left');  %,'EdgeColor','k','Margin',1
      text(0.04,0.85,stas(ista, :),'FontSize',10,'unit','normalized',...
        'horizontalalignment','left');
      print(f.fig,'-dpdf',strcat(figpath,'/aaaa_fc_gauf_zoomp1',num2str(iets),...
        strtrim(stas(ista, :)),FLAG,bndflag,'.pdf'));

      yran = [3.8 5.3];
      [f,fitobj,~,~,fgmaxfit] = plt_spgram_fcfit(tets,fgmax(:,2),fran(:,3:4),yran);
      text(0.04,0.95,num2str(years(iets)),'FontSize',10,'unit','normalized',...
        'horizontalalignment','left');  %,'EdgeColor','k','Margin',1
      text(0.04,0.85,stas(ista, :),'FontSize',10,'unit','normalized',...
        'horizontalalignment','left');
      print(f.fig,'-dpdf',strcat(figpath,'/aaaa_fc_gauf_zoomp2',num2str(iets),...
        strtrim(stas(ista, :)),FLAG,bndflag,'.pdf'));

%       ds = abs(slope1-slope2)/slope1*100;   % difference of slope in precentage
    end

  end
end   % all ets

% keyboard

%% randomly choose windows from dates after the 'active' bursts
%%% NOTE:
%%% -- 'noidates' is the dates that I consider as the totally inactive dates for the whole regions, so
%%% the waveform at the stations on these dates may result mainly from noise. These may be considered
%%% as the inter-ETS times as well (???) Difference in spectra from active tremor dates is expected. 
%%% -- 'elwhdates' is the dates that I consider as that when tremor is still active, but not mainly
%%% at the region of interest, so the waveform at the stations on these dates may result mainly from
%%% tremor sources at elsewhere. Since the main tremor front is moving generally towards strike
%%% (i.e., S -> N or SE -> NW), dates right after 'dates' should focus on the regions more northward
%%% than the region of interest, while earlier dates focus on southward regions.
noidates = [2003067;
            2003068;
            2003069;
            2004183;
            2004206;
            2004207;
            2004208;
            2005271;
            2005272;
            2005273;
            ];
          
elwhdates = [2003064;
             2003065;
             2003066;
             2004200;
             2004201;
             2004202;
             2004203;
             2005257;
             2005258;
             2005259;
             ];

focusdates =  noidates;

nfocdates = size(focusdates,1);

focwlensec = (median(trange(:,3)-trange(:,2)));
nwin = floor(86400/focwlensec);
wins = 0: focwlensec: 86400;

years = unique(floor(focusdates/1000));
nets = length(years);

% keyboard
for iets = 1: nets
  % dates in each ets
  year = years(iets);
  focdatesets = focusdates(floor(focusdates/1000)==year);
  datesets = dates(floor(dates/1000)==year);

  if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
    stas(4,:)='KELB ';
  else
    stas(4,:)='KLNB ';  % remember to change it back
  end
  
  k = 0;  % for counting
  foptcatdat = [];   % concataneted array for all windows in the same ETS, optimal component
  fortcatdat = [];   % orthogonal component
  foptsnrcatdat = [];  % optimal component within a freq range of a high SNR  
  for i = 1: length(focdatesets)
    
    %%% read daily tremor data
    date = focdatesets(i);
    jday = floor(date-year*1000);
    [STAopt,STAort,~,fileflag] = rd_daily_bpdata(year,jday,datapath,stas,PERMSTA,POLSTA,...
                              PERMROTS,POLROTS,sps,lo,hi,npo,npa,[],[],[],[]);
    
    %%% We do 2 things here:
    %%% 1. simply get the spectra of all windows and store them
    %%% 2. concatenate all windows from the same ETS as if they are continuous, then get the
    %%% spectrogram of it
    % get the windows (time ranges) of the same day, for efficiency
    rng('default');
    selwini = randi(nwin, round(size(trange,1)/nfocdates), 1);
    selwini = sort(selwini);

%     rangetemp = trange(trange(:,1)==datesets(i), :);
%     for j = 1: size(rangetemp,1)    
%       tst = rangetemp(j,2);
%       ted = rangetemp(j,3);
    for j = 1: size(selwini,1)
      tst = wins(selwini(j));
      ted = wins(selwini(j)+1);
      wlen = 2048;
      wlensec = wlen/sps;
      if ted-tst>=wlensec
        polap = 50;
        nfft = wlen;
        Fs = sps;
        foptdat = STAopt(max(floor(tst*sps),1): min(floor(ted*sps-1),86400*sps), :);
        fortdat = STAort(max(floor(tst*sps),1): min(floor(ted*sps-1),86400*sps), :);
%         foptsnrdat = STAoptsnr(max(floor(tst*sps),1): min(floor(ted*sps-1),86400*sps), :);
        for ista = 1: nsta
          [fpcopt(:,ista),fpcft]=pchave(foptdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
%           [fpcort(:,ista),~]=pchave(fortdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
%           [fpcoptsnr(:,ista),~]=pchave(foptsnrdat(:,ista+1),wlen,polap,nfft,Fs,'MAD','dpss');
        end
        
        % stores for every rtm
        k = k+1;
        fpcetsopt(:,:,k) = fpcopt(:,:);   % ifreq, ista, irtm
%         fpcetsort(:,:,k) = fpcort(:,:);
%         fpcetsoptsnr(:,:,k) = fpcoptsnr(:,:);
        
        % concatenate all windows in the same ETS
        foptcatdat = [foptcatdat; foptdat];
%         fortcatdat = [fortcatdat; fortdat];
%         foptsnrcatdat = [foptsnrcatdat; foptsnrdat];
        
        %       keyboard
      end  % if length is enough

    end % bursts in each date
    
  end % dates in each ets
  
  % store spectra of all windows of the different ETSs separately
  fpcall{iets,1} = fpcetsopt(:,:,:);
%   fpcall{iets,2} = fpcetsort(:,:,:);
%   fpcall{iets,3} = fpcetsoptsnr(:,:,:);
  
  % now for the concatenated window, get the spectrogram as if it is continuous
  %%% Recommended this way %%%%%%%%%%%%
  %%%%%%%%%%%%% the following is to test the function written by FJS
  %%% Well, it seems that 'spectrogram2' basically gives the similar result
  %%% with the same parameter setting as the built-in 'spectrogram'.
  seglen = wlen;
  %         nfft = pow2(nextpow2(seglen));
  nfft = wlen;
  window = hann(nfft);
  Fs = sps;
  olap = polap/100;
  nolap = seglen*olap;
  
  fspgcatopt = [];
  fspgcatort = [];
  fspgcatoptsnr = [];
  for ista = 1: nsta
    [fspgcatopt(:,:,ista), fspgft, foct, fpsddB]=spectrogram2(foptcatdat(:,ista+1),nfft,Fs,seglen,nolap,'s');
%     [fspgcatort(:,:,ista), ~,~,~]=spectrogram2(fortcatdat(:,ista+1),nfft,Fs,seglen,nolap,'s');
%     [fspgcatoptsnr(:,:,ista), ~,~,~]=spectrogram2(foptsnrcatdat(:,ista+1),nfft,Fs,seglen,nolap,'s');
  end
  
  % store the spectrograms of the different ETSs separately
  fspgall{iets,1} = fspgcatopt(:,:,:);
%   fspgall{iets,2} = fspgcatort(:,:,:);
%   fspgall{iets,3} = fspgcatoptsnr(:,:,:);
  ftall{iets,1} = foct;
  
end % all ets

% clean some temp variables to make more room
clear STAopt STAort 
clear foptdat fortdat 
clear fpcopt fpcort 
clear fpcetsopt fpcetsort 
clear fspgcatopt fspgcatort 

% keyboard


%% plot the spectra of optimal components of all ETS, windows from dates after the 'active' bursts
xran = [0.1 20];
yran = [1e-2 1e2];
minfnorm = 1.5; maxfnorm = 5;
[f] = plt_spectra_of_bursts_norm(years,stas,fpcft,fpcall,minfnorm,maxfnorm,xran,yran);

if isequal(focusdates,noidates)
  print(f.fig,'-dpdf',strcat(figpath,'/inactivewin_spa',FLAG,'.pdf'));
elseif isequal(focusdates,elwhdates)
  print(f.fig,'-dpdf',strcat(figpath,'/elsewhere_spa',FLAG,'.pdf'));
end
% keyboard

%% plot the spectrogram of concatanated windows for all ETS, windows from dates after the 'active' bursts 
for iets = 1: nets
  
  yran = [2 6];
  cran = [-0.4,0.4];
  minfnorm = 1.5; maxfnorm = 5;
  filtflag = 0;
  filtsigma = [];
  [f] = plt_spectrogram_of_bursts_norm(years,stas,iets,fspgft,fspgall,ftall,minfnorm,maxfnorm,yran,cran,...
                                              filtflag,filtsigma);  
  if isequal(focusdates,noidates)
    print(f.fig,'-dpdf',strcat(figpath,'/inactivewin_spg',num2str(iets),FLAG,'.pdf'));
  elseif isequal(focusdates,elwhdates)
    print(f.fig,'-dpdf',strcat(figpath,'/elsewhere_spg',num2str(iets),FLAG,'.pdf'));
  end

end   % all ets

keyboard

%% cross-correlation between spectrograms of tremor dates and other dates
%%%% this came across because the band-like feature in the spectrogram seems intriguing to us, we
%%%% are not sure about the source of contribution. It can be source or source-path, or the station.
%%%% Although it would be less interesting if it is about site effect, there is chance that some
%%%% rain or rain changes the characteristics of the site. Check the note 2021-10-20 for details
%%%% One way to think about it is, if the band feature is not seeing at the same station at other
%%%% times (when tremor is active at other regions, or when tremor is inactive so it is recording
%%%% mainly the noise). The first scenario makes sure that the source/path changes but not station.
%%%% The latter may not be ideal becasue the source might be viewed as scattered noise everywhere.

%%% 2021/10/29
%%% I am still trying to figure out the best way to do the cross-correlation. 2D CC sounds like
%%% finding a misaligned piece in a 2D space, but in our case, the time is synthesized, it is not
%%% 1-1 corresponded, and i think we care the most about what is domiant freq of the band and what
%%% is the relation between the freqs of different spectrograms. So currently i want to do a 1D CC
%%% to each 1-1 column of the spectrograms.

iets = 2; 
ista = 3;

%%% spectrogram of tremor 
tmp = spgall{iets,1};
spgtmr = squeeze(tmp(:,:,ista));
spgtmr = sqrt(spgtmr);  % convert power to amplitude

% normalize the spectrogram according to the mean amplitude within a freq range, note that each
% col of the entire array is the FFT of each overlapping window, where there is no averaging
%     minfreq = 2; maxfreq = 4; %minfreq=2.;
minfreq = 1.5; maxfreq = 5;
[~,indmin] = min(abs(spgft-minfreq));
[~,indmax] = min(abs(spgft-maxfreq));
normlzer = mean(spgtmr(indmin:indmax,:),1);
%     keyboard
for k = 1: size(spgtmr, 2)  %here k is number of overlapping subwindow for FFT, larger than number of bursts
  spgtmr(:,k)=spgtmr(:,k)/normlzer(k);
end   % all overlapping subwindows for FFT



%%% spectrogram of tremor 
tmp = fspgall{iets,1};
spgelse = squeeze(tmp(:,:,ista));
spgelse = sqrt(spgelse);  % convert power to amplitude

% normalize the spectrogram according to the mean amplitude within a freq range, note that each
% col of the entire array is the FFT of each overlapping window, where there is no averaging
%     minfreq = 2; maxfreq = 4; %minfreq=2.;
minfreq = 1.5; maxfreq = 5;
[~,indmin] = min(abs(fspgft-minfreq));
[~,indmax] = min(abs(fspgft-maxfreq));
normlzer = mean(spgelse(indmin:indmax,:),1);
%     keyboard
for k = 1: size(spgelse, 2)  %here k is number of overlapping subwindow for FFT, larger than number of bursts
  spgelse(:,k)=spgelse(:,k)/normlzer(k);
end   % all overlapping subwindows for FFT

%%% cross-correlation col by col 
lagfreq = zeros(size(spgtmr,2), 1);
maxcoef = zeros(size(spgtmr,2), 1);
for i = 1: size(spgtmr,2)
  [coef, lag] = xcorr(spgtmr(:,i), spgelse(:,i), floor(size(spgtmr,1)*0.2), 'coeff');
  [maxcoef(i), idx] = max(coef);
  lagsamp = lag(idx);
  lagfreq(i) = lagsamp*(sps/nfft);
end




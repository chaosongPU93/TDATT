% function identify_counter_pgc(fam,iup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to find the counter example of 11 (12) detections where the time lag 
% from min to max is short compared to the full duration, which is identified
% in 'identify_pgc.m', from the merged catalog of all fams of interest (LF).
%
% -- Basically, we wish to find the examples where the time lag is more 'normal'
%   and occur at roughly the same location (and a short time period), assuming
%   the attenuation would not change abruptly in space and time.
% -- The way to do it is, query the detections inside a range (rectangle/circle)
%   around the target and also occur during the same date, check their time lag
%   (in the nearby region, it might be easier to find counter examples, but maybe
%   difficult to find ones that are also close in time)  
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/12
% Last modified date:   2021/07/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
clear
close all

%% default value for easy debugging
% defval('fam', '065');
defval('iup',1);


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
temppath = strcat(datapath, '/templates/PGCtrio/');
rstpath = strcat(datapath, '/PGCtrio');
figpath = strcat(workpath, '/project2021/PGCtrio');
hypopath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forproj21');


% combine all fams into an pool, 12 fams in total, final version        
fampool = [
           '002'
           '243'; 
           '253';
           '036';
           '034'; 
           '061';
           '023';
           '251';
           '240';
           '255'; 
           '012'; 
           '065'; 
           ];    % family number
nfam = size(fampool, 1);

FLAG = 'PGC';

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];

stas=['PGC '
      'SSIB'
      'SILB'];     % determine the trio and order
nsta=size(stas,1);         %  number of stations


% load the hard copy of the high-quality detections where the time lag is pretty short
rstgood = [
    
 2004197,9402,-5,0,0.593000000000000,-0.0700000000000000,1.42800000000000,9395.50500000000,...
    0.670000000000000,0.0704200000000000,7.19231947707078,9.54071264056361,6.81044487427466,...
    11.7739508443404,5.63811048839071,14.1009211053264,6.01880341880342,9.91133004926108,0,0,0,0,240;
 2004197,45978,-11,6,0.665000000000000,-0.800000000000000,4.28600000000000,45972.5270000000,...
    0.504000000000000,0.219000000000000,24.1508601676224,5.55837563451777,25.8895850573354,...
    6.26968222158603,16.2222222222222,9.27966101694915,10.6104651162791,7.71398379711166,1,1,1,1,243;
 2004198,15858,1,5,0.510000000000000,1.60000000000000,2.10500000000000,15851.2090000000,...
    0.639000000000000,0.107900000000000,24.1333035115187,4.43850267379679,13.3358052156717,...
    4.67706978760295,9.00667779632721,6.97028423772610,9.89000916590284,9.15182357930449,0,0,0,0,243;
 2004198,84410,-17,-19,0.580000000000000,1.98000000000000,1.02800000000000,84415.1490000000,...
    0.540000000000000,0.0576400000000000,13.1418148654811,4.99480069324090,14.1274509803922,...
    4.18895348837209,14.9287749287749,5.43261074458058,11.1792086889061,4.19505094614265,0,0,0,0,240;
 2004199,906,1,4,0.710000000000000,-0.610000000000000,12,900.351000000000,0.803000000000000,...
    0.672300000000000,55.0614250614251,19.1266002844950,68.9821465216499,23.3275503122831,...
    50.4351087771943,25.5336118496012,23.6975678533662,35.1437532671197,0,1,0,1,243;
 2004199,2306,-8,1,0.686000000000000,0.840000000000000,1.82300000000000,2305.50100000000,...
    0.537000000000000,0.0982200000000000,33.0929919137466,9.37213740458015,20.3818219547624,...
    11.1121167552891,10.6981810260320,11.6789536266350,5.93115942028986,11.5620953502060,0,1,0,0,243;
 2005255,19242,-12,4,0.672000000000000,0.810000000000000,7.32600000000000,19246.5610000000,...
    0.689000000000000,0.375000000000000,22.6176115802171,8.85896527285613,23.7793278376665,...
    8.79043600562588,14.6255850234009,10.7051099057950,17.9856115107914,12.8556736372986,0,1,0,1,243;
 2005255,41338,14,0,0.727000000000000,0.0100000000000000,2.42700000000000,41344.9740000000,...
    0.767000000000000,0.132500000000000,26.8817204301075,9.79305247597931,28.2696820994239,...
    11.2862010221465,23.2741963815212,15.2597028676725,24.1479861490796,14.2167381974249,0,0,0,1,2;
 2005255,67938,0,-8,0.657000000000000,-0.650000000000000,2.49700000000000,67934.4890000000,...
    0.549000000000000,0.127100000000000,15.2489502099580,5.98117647058824,12.8552644887226,...
    7.12044817927171,7.80712530712531,9.01418439716312,8.23720025923526,10.3670473083197,0,1,1,1,2;
 2005255,72234,1,5,0.714000000000000,0.450000000000000,2.79400000000000,72240.9920000000,...
    0.831000000000000,0.145400000000000,11.8211382113821,9.52193844138834,14.5400000000000,...
    11.6413130504404,23.8634498604957,17.6263789550249,21.4580873671783,21.0572049239681,0,1,1,1,243;
 2005256,13354,3,2,0.789000000000000,-0.530000000000000,2.89100000000000,13352.7080000000,...
    0.615000000000000,0.166800000000000,22.6538095884830,7.49326145552561,21.3081246806336,...
    9.79448032883147,30.9691793538804,17.1516709511568,16.3209393346380,22.7433869648214,0,1,1,1,243;
 ];


%%% load locations of isolated lf detections
fname = strcat(hypopath, '/evtloc.offset_all_isolatedlf');
isolf = load(fname);
[dx, dy] = absloc2relaloc(isolf(:,1),isolf(:,2),-123.585000, 48.436667);
isolf = [dx dy isolf];

isolfsort = [isolf(11,:);
             isolf(3:4,:);
             isolf(12,:);
             isolf(5:7,:);
             isolf(1:2,:);
             isolf(8:9,:);
             ];
             
% paste the information
rstgood = [isolfsort rstgood];


%%% load the merged catalog of hf and lf detections 
winlensechf = 4;
winlensec = 16;
loopoffmax = 4;
xcmaxAVEnmin = 0.45;

SUFFIXhf = strcat('up.hf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
                  '.',num2str(xcmaxAVEnmin));
SUFFIXlf = strcat('lf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
                  '.',num2str(xcmaxAVEnmin));
hffname = strcat(hypopath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXhf);
lffname = strcat(hypopath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXlf);

hftime = load(hffname);
lftime = load(lffname);

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


%% 
for i = 9: size(rstgood, 1)
    
    % boundary of a circular region around the target high-quality detections
    x0 = rstgood(i, 1);
    y0 = rstgood(i, 2);
    rad = 2.5;
    ddeg = 10;
    [x, y] = circle_chao(x0,y0,rad,ddeg);
    bnd = [x' y';
           x(1) y(1);
           ];
    
    % events inside bound
    [iin,ion] = inpolygon(lftime(:,1),lftime(:,2),bnd(:,1),bnd(:,2));
    isinbnd = iin | ion;
    lfbnd = lftime(isinbnd == 1, :);
    
    % target date
    date = rstgood(i, 10);
    
    % events inside bound and also occur in the same date
    lfbndtime = lfbnd(lfbnd(:,14) == date, :);
    
    %%% for these eligible detections, let's see if we can find any counter example
    if isempty(lfbndtime)
        sprintf('No detections close in space and time with the target %d',i);
    else
        sps = 20;
        winlensec=16;       % CHECK in detection script for what were used in first-year report
        winoffsec=8;        % window offset in sec, which is the step of a moving window
        winlen=winlensec*sps;      % length in smaples
        hi=1.25;
        lo=0.5;
        npo=2;     % poles, passes of filters
        npa=2;
        concentration=1.1; %in seconds; how concentrated is the coherent energy within the window?
        cncntr=concentration*sps;   % in samples, 20
        
        stasnew=['LZB  '
            'TWKB '
            'MGCB '
            'KLNB '];  % twkb lzb mgcb
        nstasnew=size(stasnew,1);
        
        for j = 1: size(lfbndtime,1)
            famnum = lfbndtime(j, 13);
            if famnum <= 9
                fam=['00',int2str(famnum)];
            elseif famnum<= 99
                fam=['0',int2str(famnum)];
            else
                fam=int2str(famnum);
            end
            year = floor(date/1000);
            jday = floor(date-year*1000);
            dates = jul2dat(year,jday);
            YEAR = int2str(dates(3));
            DAY = int2str(dates(2));
            if jday <= 9
                JDAY=['00',int2str(jday)];
            elseif jday<= 99
                JDAY=['0',int2str(jday)];
            else
                JDAY=int2str(jday);
            end
            if dates(1) == 7
                MO = {'Jul. '};
            elseif dates(1) == 3
                MO = {'Mar. '};
            elseif dates(1) == 9
                MO = {'Sep. '};
            end
            
            if isequal(fam,'002')
                cyclskip = 20;
                mshift=19+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
                loopoffmax=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
                xcmaxAVEnmin=0.5; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
            else
                cyclskip = 20;
                mshift=19+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
                loopoffmax=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
                xcmaxAVEnmin=0.45; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
            end
            
            IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),...
                '.',int2str(npo),int2str(npa),'.ms', int2str(mshift)];
            
            % maybe consider to load the complete 1-day trace and truncate a segment of the same length and
            % is centered at the main arrival, the advantage is, visually every detection is at the center,
            % and 'findpeaks' would run without any issue when finding the nearest local max and min. But it
            % can also misleads people and myself to think that is the detecting window, which is not true,
            % as the detection might be located in anywhere within the real detecting window (even if i
            % saved the duplicates that are closest to the center of the window
            fname2 = strcat(rstpath, '/MAPS/pj21trace1d_',IDENTIF,'_',num2str(lo),'-',num2str(hi),...
                '_',num2str(winlen/sps),'s',num2str(sps),'sps');
            trace1df = load(fname2);
            % trace1df = [timsSTA' STAopt'];
            
            % load the no filtered trace of one day, this would help solve the issue that record length is
            % shorter than that required by the longest period of the filter
            % Note that there is no pre-alignment in this case
            fname = strcat(rstpath, '/MAPS/pj21trace1dbb_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                num2str(winlen/sps),'s',num2str(sps),'sps');
            trace1dbbf = load(fname);
            
            %%% load the template filtering effect result
            lolf = 0.5;
            hilf = 1.25;
            lohf = 1.25;
            hihf = 6.5;
            winsechf = 4;
            winseclf = winlensec;
            if isequal(fam,'002')
                lofflf = 4;
                ccminlf = 0.5;
            else
                lofflf = 4;
                ccminlf = 0.45;
            end
            PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
                '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(lofflf),'.ccm', ...
                num2str(ccminlf),'.','80sps');
            fname = strcat(rstpath, '/MAPS/tempfeff_PGC_enc','_',PREFIX);
            filtcor = load(fname);
            spsratio = sps/80;
            filtlf = filtcor(:,2)*spsratio;   % hf template shift due to filterig effect, sign is the same

            %%%%%%%%%%%%%%%%% PLOT 2, octave-filtered version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f.fig = figure;
            f.fig.Renderer = 'painters';
            widin = 6;  % maximum width allowed is 8.5 inches
            htin = 9.5;   % maximum height allowed is 11 inches
            set(f.fig,'Position',[1*scrsz(3)/8 scrsz(4)/10 widin*res htin*res]);
            nrow = 9;
            ncol = 1;
            for isub = 1:nrow*ncol
                f.ax(isub) = subplot(nrow,ncol,isub);
                f.ax(isub).FontSize = 8;
                f.ax(isub).Box = 'on';
                grid(f.ax(isub), 'on');
                f.ax(isub).GridLineStyle = '--';
            end
            
            % reposition
            %         set(f.ax(11), 'position', [ 0.1, 0.06, 0.85, 0.08]);
            %         set(f.ax(10), 'position', [ 0.1, 0.14, 0.85, 0.08]);
            %         set(f.ax(9), 'position', [ 0.1, 0.22, 0.85, 0.08]);
            set(f.ax(9), 'position', [ 0.1, 0.05, 0.85, 0.21]);     % for spectrum
            %         set(f.ax(9), 'position', [ 0.1, 0.05, 0.34, 0.21]);      % for amp variation
            set(f.ax(8), 'position', [ 0.1, 0.30, 0.85, 0.08]);     % 8-4, for octave filtering
            set(f.ax(7), 'position', [ 0.1, 0.38, 0.85, 0.08]);
            set(f.ax(6), 'position', [ 0.1, 0.46, 0.85, 0.08]);
            set(f.ax(5), 'position', [ 0.1, 0.54, 0.85, 0.08]);
            set(f.ax(4), 'position', [ 0.1, 0.62, 0.85, 0.08]);
            set(f.ax(3), 'position', [ 0.1, 0.70, 0.85, 0.08]);     % for broader band
            set(f.ax(2), 'position', [ 0.1, 0.78, 0.85, 0.08]);     % for HF
            set(f.ax(1), 'position', [ 0.1, 0.86, 0.85, 0.12]);     % for LF
            
            %         WAVE = 'win';
            WAVE = 'day';
            %         TAPER = 'cosine';
            %         TAPER = 'hann';
            TAPER = 'tukey';
            %         TAPER = 'none';
            
            if strcmp(WAVE, 'day')
                %%%%% if use one-day trace rather than the detecting window
                minfreq = -2;   % lowest oxtave frequency power, ie, lowest freq is 2^-(minfreq)
                nwlensec = max(4*pow2(-minfreq), 32);
                nwinlen = nwlensec*sps;
                %             nwinlen = 64*sps;
                [~,ind] = min(abs(lfbndtime(j, 16)-trace1df(:,1)));
                % here, instead of using 'idiff', we use idiff+ cncntr/2 as a rough estimate to 0-crossing
                ind = ind+cncntr/2;
                trace1d = [];
                trace1d(:,1:2) = trace1df(ind-nwinlen/2: ind+nwinlen/2-1, 1:2);
                % although no pre-alignment, but aligned here according to the results
                off12fec = lfbndtime(j, 9);
                off13fec = lfbndtime(j, 10);
                filt12 = filtlf(1);
                filt13 = filtlf(2);
                off12 = off12fec + filt12;
                off13 = off13fec + filt13;
                trace1d(:,3) = trace1df(ind-nwinlen/2 - off12: ind+nwinlen/2-1 - off12, 3);
                trace1d(:,4) = trace1df(ind-nwinlen/2 - off13: ind+nwinlen/2-1 - off13, 4);
                trace = trace1d;
                
                [~,ind] = min(abs(lfbndtime(j, 16)-trace1dbbf(:,1)));
                % here, instead of using 'idiff', we use idiff+ cncntr/2 as a rough estimate to 0-crossing
                ind = ind+cncntr/2;
                trace1dbb = [];
                % +/-1 sec for alignment in different passbands
                trace1dbb(:,1:2) = trace1dbbf(ind-nwinlen/2-sps: ind+nwinlen/2+sps-1, 1:2);
                trace1dbb(:,3) = trace1dbbf(ind-nwinlen/2-sps - off12: ind+nwinlen/2+sps-1 - off12, 3);
                trace1dbb(:,4) = trace1dbbf(ind-nwinlen/2-sps - off13: ind+nwinlen/2+sps-1 - off13, 4);
                tracebb = trace1dbb;
                
                idiff = nwinlen/2-cncntr/2;
                %%%%%
            end
            
            tracebbtap = tracebb;
            if strcmp(TAPER, 'cosine')
                % cosine taper over 1 second before filtering:
                aaa = 0:pi/(2*sps):pi/2-pi/(2*sps);
                w = sin(aaa');
                tracebbtap(1:sps,2:4) = w.*tracebb(1:sps,2:4);
                tracebbtap(end-sps+1:end,2:4) = flipud(w).*tracebb(end-sps+1:end,2:4);
            elseif strcmp(TAPER, 'hann')
                % hann taper, in fact just a hann window
                w = hann(nwinlen);
                tracebbtap(:,2:4) = w.* tracebb(:,2:4);
            elseif strcmp(TAPER, 'none')
                tracebbtap = tracebb;
            end
            
            %%% LF
            ax = f.ax(1);
            hold(ax,'on');
            
            [maxamppos,indmax] = max(trace(max(nwinlen/2-cncntr,1):min(nwinlen/2+cncntr-1,nwinlen),...
                2:4), [],1);
            [minampneg,indmin] = min(trace(max(nwinlen/2-cncntr,1):min(nwinlen/2+cncntr-1,nwinlen),...
                2:4), [],1);
            
            %roughly estimate the time of zero-crossing which indicate the origin time of a signal.
            indmax = indmax-1+max(nwinlen/2-cncntr,1);  % convert to global index
            indmin = indmin-1+max(nwinlen/2-cncntr,1);
            %in case indmin is not smaller than indmax for some detections
            if indmin(1) < indmax(1)
                seg = trace(indmin(1): indmax(1), 2);  % for zero-crossing timing, only use the main station
                [~,zc] = min(abs(seg));
                zc = zc-1+indmin(1);  % convert to global index
            else
                seg = trace(indmax(1): indmin(1), 2);  % for zero-crossing timing, only use the main station
                [~,zc] = min(abs(seg));
                zc = zc-1+indmax(1);  % convert to global index
            end
            is = trace(1, 1);   % start time of each win
            ien= trace(end, 1);         % end time of each win
            zctime = (is+zc/sps);
            
            ymami = max([maxamppos; -minampneg], [], 1);
            pltscale = 1.2;
            ym = pltscale*max(ymami);
            axis(ax,[is ien -ym ym]);
            
            %%% patch the main arrival window
            xarea=[is+idiff/sps is+(idiff+cncntr-1)/sps is+(idiff+cncntr-1)/sps is+idiff/sps ...
                is+idiff/sps];
            yarea=[-ym -ym ym ym -ym];
            patch(ax,xarea,yarea,'k','Facealpha',0.15,'edgecolor','none');
            
            % annotate the start and end of the original detecting window
            ocnt = lfbndtime(j, 15);
            ois = ocnt-winlensec/2;
            oien = ocnt+winlensec/2;
            plot(ax,[ois ois],ax.YLim,'k:','linew',1.5);
            plot(ax,[oien oien],ax.YLim ,'k:','linew',1.5);
            
            %%% plot the seismogram at 3 stations of each win
            plot(ax,trace(:, 1),trace(:, 2),'r','linew',0.5);
            plot(ax,trace(:, 1),trace(:, 3),'b','linew',0.5);
            plot(ax,trace(:, 1),trace(:, 4),'k','linew',0.5);
            
            %%% plot freq. band
            text(ax,0.02,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
            text(ax,0.02,0.8,strcat({'Aligned'}),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
            text(ax,0.98,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','r');
            text(ax,0.98,0.8,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','b');
            text(ax,0.98,0.7,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','k');
            text(ax,0.02,0.1,sprintf('%s  %s',FLAG,fam),'fontsize',9,'unit','normalized',...
                'Horizontalalignment','left','EdgeColor','k','Margin',0.5);
            text(ax,0.35, 0.9, num2str(off13),'fontsize',8,'color','k','unit','normalized',...
                'Horizontalalignment','left');
            text(ax,0.25, 0.9, num2str(off12),'fontsize',8,'color','b','unit','normalized',...
                'Horizontalalignment','left');
            text(ax,0.75,0.9,sprintf('%.4f',lfbndtime(j, 17)),'fontsize',8,'unit','normalized',...
                'Horizontalalignment','right');
            % if it is also confirmed to be real by the 4th station
            %         flagconf = [];
            for jj = 1: nstasnew
                if lfbndtime(j, 43+jj-1) == 1
                    flagconf(jj) = strcat(stasnew(jj,1));
                else
                    flagconf(jj) = ' ';
                end
            end
            %         if allrst_new(ndet, 22)==1 || (sum(allrst_new(ndet, 20:21))>=1 && allrst_new(ndet, 22)==0)
            %             text(ax,0.98,0.1,'Confirmed','fontsize',9,'unit','normalized',...
            %                  'Horizontalalignment','right','EdgeColor','k','Margin',0.5);
            %         end
            text(ax,0.98,0.1,strcat(flagconf(1),flagconf(2),flagconf(3),flagconf(4)),'fontsize',9,...
                'unit','normalized','Horizontalalignment','right','EdgeColor','k','Margin',0.5);
            text(ax,0.25,0.1,sprintf('%.2f s',zctime),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
            
            set(ax,'XTick',is:4:ien,'fontsize',8);
            ax.XTickLabel = [];
            ylabel(ax,'Amplitude','fontsize',9);
            %         xlabel(ax,'Time (s)','fontsize',9);
            [waveenvup, waveenvlo] = envelope(trace(:,2:4));
            medenvup = median(waveenvup,2);
            medenvlo = median(waveenvlo,2);
            plot(ax,trace(:, 1),medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
            plot(ax,trace(:, 1),medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
            %         plot(ax, ax.XLim, [0.5*maxhil 0.5*maxhil], '--');
            %         plot(ax, ax.XLim, [0.25*maxhil 0.25*maxhil], '--');
            %         plot(ax, ax.XLim, [0.125*maxhil 0.125*maxhil], '--');
            hold(ax,'off');
            
            %%%%%% original HF
            ax = f.ax(2);
            hold(ax,'on');
            
            lohf = 1.25;
            hihf = 6.5;
            npo = 2;
            npa = 2;
            
            if strcmp(TAPER, 'tukey')
                % taper with tukeywin, which is actually a tapered cosine window
                %tapered length is adaptative to frequency, maybe at least longer than one full period length of
                %the lowest frequency
                fractap = round(1/lohf*sps)/size(tracebb,1)*2;
                fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
                w = tukeywin(size(tracebb,1),fractap);
                tracebbtap(:,2:4) = w.* tracebb(:,2:4);
            end
            
            tracebpadd = tracebbtap;
            for ista = 1: nsta
                tracebpadd(:,ista+1) = Bandpass(tracebbtap(:,ista+1), sps, lohf, hihf, npo, npa, 'butter');
            end
            
            % re-align the filtered traces, because the alignment was in LF, which may not be true for
            % higher frequency
            mshiftadd = 5;    % maximum allowed shift between 2 traces
            mid = ceil(size(tracebpadd,1)/2);
            fixlen = max(2*sps, 1/lohf*sps);
            loffmax = 4;
            ccmin = 0.3;  % 0.4/0.35
            iup = 1;    % times of upsampling
            BETALIGN = 1;
            %         [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc(tracebpadd(:,2:4)',mid,fixlen,...
            %                                                                    mshift,loffmax,ccmin);
            [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd(:,2:4)',mid,...
                fixlen,mshiftadd,loffmax,ccmin,iup);
            % if a better alignment cannot be achieved, use the LF alignment
            if off12add == mshiftadd+1 && off13add == mshiftadd+1
                off12add = 0;
                off13add = 0;
                BETALIGN = 0;
            end
            
            % align the records
            istart = sps+1;
            iend = nwinlen+sps;
            tracebp = [];
            tracebp(:, 1:2) = tracebpadd(istart: iend, 1:2);
            tracebp(:, 3) = tracebpadd(istart-round(off12add): iend-round(off12add), 3);
            tracebp(:, 4) = tracebpadd(istart-round(off13add): iend-round(off13add), 4);
            
            [maxamppos,indmax] = max(tracebp(max(nwinlen/2-0.5*fixlen,1): ...
                min(nwinlen/2+0.5*fixlen-1,nwinlen),2:4), [],1);
            [minampneg,indmin] = min(tracebp(max(nwinlen/2-0.5*fixlen,1): ...
                min(nwinlen/2+0.5*fixlen-1,nwinlen),2:4), [],1);
            if BETALIGN
                %roughly estimate the time of zero-crossing which indicate the origin time of a signal.
                indmax = indmax-1+max(nwinlen/2-0.5*fixlen,1);  % convert to global index
                indmin = indmin-1+max(nwinlen/2-0.5*fixlen,1);
                %in case indmin is not smaller than indmax for some detections
                if indmin(1) < indmax(1)
                    seg = tracebp(indmin(1): indmax(1), 2);  % for zero-crossing timing, only use the main station
                    [~,zc] = min(abs(seg));
                    zc = zc-1+indmin(1);  % convert to global index
                else
                    seg = tracebp(indmax(1): indmin(1), 2);  % for zero-crossing timing, only use the main station
                    [~,zc] = min(abs(seg));
                    zc = zc-1+indmax(1);  % convert to global index
                end
                zctime = (is+zc/sps);
            end
            
            ymami = max([maxamppos; -minampneg], [], 1);
            ym = pltscale*max(ymami);
            axis(ax,[is ien -ym ym]);
            
            %annotate the window to find the max/min
            plot(ax,[is+(nwinlen/2-0.5*fixlen)/sps is+(nwinlen/2-0.5*fixlen)/sps],ax.YLim,'k--',...
                'linew',1.5);
            plot(ax,[is+(nwinlen/2+0.5*fixlen-1)/sps is+(nwinlen/2+0.5*fixlen-1)/sps],ax.YLim,'k--',...
                'linew',1.5);
            
            %%% plot the seismogram at 3 stations of each win
            plot(ax,tracebp(:, 1),tracebp(:, 2),'r','linew',0.5);
            plot(ax,tracebp(:, 1),tracebp(:, 3),'b','linew',0.5);
            plot(ax,tracebp(:, 1),tracebp(:, 4),'k','linew',0.5);
            
            %%% plot freq. band
            text(ax,0.02,0.9,strcat(num2str(lohf),'-',num2str(hihf),{' Hz'}),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
            if BETALIGN
                text(ax,0.02,0.75,strcat({'Aligned'}),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','left');
                text(ax,0.25,0.1,sprintf('%.2f s',zctime),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','left');
            else
                text(ax,0.02,0.75,strcat({'Aligned (poor)'}),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','left');
            end
            text(ax,0.98,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','r');
            text(ax,0.98,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','b');
            text(ax,0.98,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','k');
            text(ax,0.35, 0.9, num2str(off13+off13add),'fontsize',8,'color','k','unit','normalized',...
                'Horizontalalignment','left');
            text(ax,0.25, 0.9, num2str(off12+off12add),'fontsize',8,'color','b','unit','normalized',...
                'Horizontalalignment','left');
            text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
                'right');
            
            set(ax,'XTick',is:4:ien,'fontsize',8);
            ax.XTickLabel = [];
            %         ylabel(ax,'Amplitude','fontsize',9);
            %         ax.YTick = [];
            hold(ax,'off');
            
            
            %%%%%%%% broader band, instead of broadband
            ax = f.ax(3);
            hold(ax,'on');
            
            lohf = 0.5;
            hihf = 6.5;
            npo = 2;
            npa = 2;
            
            if strcmp(TAPER, 'tukey')
                % taper with tukeywin, which is actually a tapered cosine window
                %tapered length is adaptative to frequency, maybe at least longer than one full period length of
                %the lowest frequency
                fractap = round(1/lohf*sps)/size(tracebb,1)*2;
                fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
                w = tukeywin(size(tracebb,1),fractap);
                tracebbtap(:,2:4) = w.* tracebb(:,2:4);
            end
            
            tracebpadd = tracebbtap;
            for ista = 1: nsta
                tracebpadd(:,ista+1) = Bandpass(tracebbtap(:,ista+1), sps, lohf, hihf, npo, npa, 'butter');
            end
            
            % re-align the filtered traces, because the alignment was in LF, which may not be true for
            % higher frequency
            %         mshift = 10;    % maximum allowed shift between 2 traces
            mid = ceil(size(tracebpadd,1)/2);
            fixlen = max(2*sps, 1/lohf*sps);
            %         loffmax = 4;
            %         ccmin = 0.3;  % 0.4/0.35
            %         iup = 1;    % times of upsampling
            BETALIGN = 1;
            [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd(:,2:4)',mid,...
                fixlen,mshiftadd,loffmax,ccmin,iup);
            
            % if a better alignment cannot be achieved, use the LF alignment
            if off12add == mshiftadd+1 && off13add == mshiftadd+1
                off12add = 0;
                off13add = 0;
                BETALIGN = 0;
            end
            
            % align the records
            istart = sps+1;
            iend = nwinlen+sps;
            tracebp = [];
            tracebp(:, 1:2) = tracebpadd(istart: iend, 1:2);
            tracebp(:, 3) = tracebpadd(istart-round(off12add): iend-round(off12add), 3);
            tracebp(:, 4) = tracebpadd(istart-round(off13add): iend-round(off13add), 4);
            
            [maxamppos,indmax] = max(tracebp(max(nwinlen/2-0.5*fixlen,1): ...
                min(nwinlen/2+0.5*fixlen-1,nwinlen),2:4), [],1);
            [minampneg,indmin] = min(tracebp(max(nwinlen/2-0.5*fixlen,1): ...
                min(nwinlen/2+0.5*fixlen-1,nwinlen),2:4), [],1);
            if BETALIGN
                %roughly estimate the time of zero-crossing which indicate the origin time of a signal.
                indmax = indmax-1+max(nwinlen/2-0.5*fixlen,1);  % convert to global index
                indmin = indmin-1+max(nwinlen/2-0.5*fixlen,1);
                %in case indmin is not smaller than indmax for some detections
                if indmin(1) < indmax(1)
                    seg = tracebp(indmin(1): indmax(1), 2);  % for zero-crossing timing, only use the main station
                    [~,zc] = min(abs(seg));
                    zerocrs1 = find(diff(sign(seg)));
                    zc = zc-1+indmin(1);  % convert to global index
                else
                    seg = tracebp(indmax(1): indmin(1), 2);  % for zero-crossing timing, only use the main station
                    [~,zc] = min(abs(seg));
                    zc = zc-1+indmax(1);  % convert to global index
                end
                zctime = (is+zc/sps);
                
                %%%compare the time lag from the min to max, with the duration of the dipole
                tlag = zeros(1,3);
                zclt = zeros(1,3);
                zcrt = zeros(1,3);
                tdura = zeros(1,3);
                for ii = 1: size(indmax,2)
                    %in case indmin is not smaller than indmax for some detections
                    if indmin(ii) < indmax(ii)
                        %roughly estimate the time lag from the min to max
                        tlag(ii) = (indmax(ii)-indmin(ii))/sps;
                        %what is the rough duration of the signal?
                        %we need to find the nearest 2 zeros on the left and right of the zero-crossing
                        seg = tracebp(indmin(ii)-sps: indmin(ii)-1, ii+1);
                        tmp = find(diff(sign(seg)),1,'last');
                        if ~isempty(tmp)
                            zclt(ii) = tmp;
                            zclt(ii) = zclt(ii)-1+indmin(ii)-sps;  % convert to global index
                        else
                            zclt(ii) = indmin(ii)-sps;
                        end
                        seg = tracebp(indmax(ii)+1: indmax(ii)+sps, ii+1);
                        tmp = find(diff(sign(seg)),1,'first');
                        if ~isempty(tmp)
                            zcrt(ii) = tmp;
                            zcrt(ii) = zcrt(ii)+indmax(ii);  % convert to global index
                        else
                            zcrt(ii) = indmax(ii)+sps;
                        end
                        tdura(ii) = (zcrt(ii)-zclt(ii))/sps;
                    else
                        tlag(ii) = (indmin(ii)-indmax(ii))/sps;
                        seg = tracebp(indmax(ii)-sps: indmax(ii)-1, ii+1);
                        tmp = find(diff(sign(seg)),1,'last');
                        if ~isempty(tmp)
                            zclt(ii) = tmp;
                            zclt(ii) = zclt(ii)-1+indmax(ii)-sps;  % convert to global index
                        else
                            zclt(ii) = indmax(ii)-sps;
                        end
                        seg = tracebp(indmin(ii)+1: indmin(ii)+sps, ii+1);
                        tmp = find(diff(sign(seg)),1,'first');
                        if ~isempty(tmp)
                            zcrt(ii) = tmp;
                            zcrt(ii) = zcrt(ii)+indmin(ii);  % convert to global index
                        else
                            zcrt(ii) = indmin(ii)+sps;
                        end
                        tdura(ii) = (zcrt(ii)-zclt(ii))/sps;
                    end
                end
                trat = tlag./tdura;
                
            end
            
            ymami = max([maxamppos; -minampneg], [], 1);
            ym = pltscale*max(ymami);
            axis(ax,[is ien -ym ym]);
            
            %annotate the window to find the max/min
            plot(ax,[is+(nwinlen/2-0.5*fixlen)/sps is+(nwinlen/2-0.5*fixlen)/sps],ax.YLim,'k--',...
                'linew',1.5);
            plot(ax,[is+(nwinlen/2+0.5*fixlen-1)/sps is+(nwinlen/2+0.5*fixlen-1)/sps],ax.YLim,'k--',...
                'linew',1.5);
            
            %%% plot the seismogram at 3 stations of each win
            plot(ax,tracebp(:, 1),tracebp(:, 2),'r','linew',0.5);
            plot(ax,tracebp(:, 1),tracebp(:, 3),'b','linew',0.5);
            plot(ax,tracebp(:, 1),tracebp(:, 4),'k','linew',0.5);
            
            %         text(ax,0.02,0.9,strcat('Broadband'),'fontsize',9,'unit',...
            %              'normalized','Horizontalalignment','left');
            text(ax,0.02,0.9,strcat(num2str(lohf),'-',num2str(hihf),{' Hz'}),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
            if BETALIGN
                text(ax,0.02,0.75,strcat({'Aligned'}),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','left');
                text(ax,0.25,0.1,sprintf('%.2f s',zctime),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','left');
                text(ax,0.75,0.1,sprintf('%.2f/%.2f=%.2f',median(tlag),median(tdura),median(trat)),...
                    'fontsize',8,'unit','normalized','Horizontalalignment','right');
                
            else
                text(ax,0.02,0.75,strcat({'Aligned (poor)'}),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','left');
            end
            text(ax,0.98,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','r');
            text(ax,0.98,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','b');
            text(ax,0.98,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','right','color','k');
            text(ax,0.35, 0.9, num2str(off13+off13add),'fontsize',8,'color','k','unit','normalized',...
                'Horizontalalignment','left');
            text(ax,0.25, 0.9, num2str(off12+off12add),'fontsize',8,'color','b','unit','normalized',...
                'Horizontalalignment','left');
            text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
                'right');
            
            set(ax,'XTick',is:4:ien,'fontsize',8);
            %         ylabel(ax,'Amplitude','fontsize',9);
            ax.XTickLabel = [];
            %         ax.YTick = [];
            hold(ax,'off');
            
            
            %%%%%%% HF, octave-passband filtered
            octavef = 2.^(3:-1:minfreq);
            ymmat = zeros(length(octavef),3); %the matrix storing the max/min of the octave-filtered main arrival
            ymmat(1, :) = ymami;    %this is from the broadband above
            efffreq = [];
            for ifreq = 1: length(octavef)-1
                
                ax = f.ax(ifreq+3);
                hold(ax,'on');
                
                lohf = octavef(ifreq+1);
                hihf = octavef(ifreq);
                npo = 2;    % npo being 2 or 3 is proper
                npa = 2;
                
                if strcmp(TAPER, 'tukey')
                    % taper with tukeywin, which is actually a tapered cosine window
                    %tapered length is adaptative to frequency, maybe at least longer than one full period length of
                    %the lowest frequency
                    fractap = round(1/lohf*sps)/size(tracebb,1)*2;
                    fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
                    w = tukeywin(size(tracebb,1),fractap);
                    tracebbtap(:,2:4) = w.* tracebb(:,2:4);
                end
                
                tracebpadd = tracebbtap;
                if diff(2*[lohf hihf]/sps)>=0.01  % the bandwidth is moderate
                    for ista = 1: nsta
                        tracebpadd(:,ista+1) = Bandpass(tracebbtap(:,ista+1),sps,lohf,hihf,npo,npa,'butter');
                    end
                else   % the bandwidth is too narrow, thus the filtering may be less reliable
                    nsps = 5;
                    [num, denom] = rat(nsps/sps);
                    tracebbtapds = zeros(size(tracebbtap,1)*num/denom, nsta);
                    tracedsbp = zeros(size(tracebbtap,1)*num/denom, nsta);
                    for ista = 1: nsta
                        %downsample first
                        tracebbtapds(:,ista) = resample(tracebbtap(:,ista+1),num,denom);   % times = num/denom
                        %filter then
                        tracedsbp(:,ista) = Bandpass(tracebbtapds(:,ista),nsps,lohf,hihf,npo,npa,'butter');
                        %interpolate finally
                        tracebpadd(:,ista+1) = interp(tracedsbp(:,ista),denom/num);
                    end
                    text(ax,0.16,0.9,strcat({'narrow passband'}),'fontsize',9,'unit',...
                        'normalized','Horizontalalignment','left');
                end
                
                % re-align the filtered traces, because the alignment was in LF, which may not be true for
                % higher frequency
                %             mshift = 10;    % maximum allowed shift between 2 traces
                mid = ceil(size(tracebpadd,1)/2);
                fixlen = max(2*sps, 1/lohf*sps);
                %             loffmax = 4;
                %             ccmin = 0.3;  % 0.4/0.35
                %             iup = 1;    % times of upsampling
                BETALIGN = 1;
                [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd(:,2:4)',mid,...
                    fixlen,mshiftadd,loffmax,ccmin,iup);
                % if a better alignment cannot be achieved, use the LF alignment
                if off12add == mshiftadd+1 && off13add == mshiftadd+1
                    off12add = 0;
                    off13add = 0;
                    BETALIGN = 0;
                end
                
                % align the records
                istart = sps+1;
                iend = nwinlen+sps;
                tracebp = [];
                tracebp(:, 1:2) = tracebpadd(istart: iend, 1:2);
                tracebp(:, 3) = tracebpadd(istart-round(off12add): iend-round(off12add), 3);
                tracebp(:, 4) = tracebpadd(istart-round(off13add): iend-round(off13add), 4);
                
                [maxamppos,indmax] = max(tracebp(max(nwinlen/2-0.5*fixlen,1): ...
                    min(nwinlen/2+0.5*fixlen-1,nwinlen),2:4), [],1);
                [minampneg,indmin] = min(tracebp(max(nwinlen/2-0.5*fixlen,1): ...
                    min(nwinlen/2+0.5*fixlen-1,nwinlen),2:4), [],1);
                if BETALIGN
                    %roughly estimate the time of zero-crossing which indicate the origin time of a signal.
                    indmax = indmax-1+max(nwinlen/2-0.5*fixlen,1);  % convert to global index
                    indmin = indmin-1+max(nwinlen/2-0.5*fixlen,1);
                    %in case indmin is not smaller than indmax for some detections
                    if indmin(1) < indmax(1)
                        seg = tracebp(indmin(1): indmax(1), 2);  % for zero-crossing timing, only use the main station
                        [~,zc] = min(abs(seg));
                        zc = zc-1+indmin(1);  % convert to global index
                    else
                        seg = tracebp(indmax(1): indmin(1), 2);  % for zero-crossing timing, only use the main station
                        [~,zc] = min(abs(seg));
                        zc = zc-1+indmax(1);  % convert to global index
                    end
                    zctime = (is+zc/sps);
                end
                
                ymami = max([maxamppos; -minampneg], [], 1);
                ym = pltscale*max(ymami);
                ymmat(ifreq+1, :) = ymami;
                
                axis(ax,[is ien -ym ym]);
            
                % the octave-filtered can show you the amplitude of main arrival vs. the surrounding
                % time (i.e., ideally noise or mixture with other arrivals), and the coherency
                % across stations, these information may tell us which range of frequencies can we
                % trust on the spectrum!
                [waveenvup, waveenvlo] = envelope(tracebp(:,2:4));
                medenvup = median(waveenvup,2);
                medenvlo = median(waveenvlo,2);
                medfull = median(mean(medenvup-medenvlo));  % median amp of the full trace
                medfix = median(ymami); % median of the max amp with the fixed len centered at main arrival
                
                % -- we need the passband to have a good alignment, reflected by 'BETALIGN' flag, or
                % the CC coef 'cc'
                % -- we also need 'medfix' to be at least higher than 'medfull'
                if BETALIGN && medfix >= medfull
                    efffreq = [efffreq lohf hihf];
                end
                efffreqhi = max(efffreq);
                efffreqlo = min(efffreq);
                
                %annotate the window to find the max/min
                plot(ax,[is+(nwinlen/2-0.5*fixlen)/sps is+(nwinlen/2-0.5*fixlen)/sps],ax.YLim,'k--',...
                    'linew',1.5);
                plot(ax,[is+(nwinlen/2+0.5*fixlen-1)/sps is+(nwinlen/2+0.5*fixlen-1)/sps],ax.YLim,'k--',...
                    'linew',1.5);
                
                %%% plot the seismogram at 3 stations of each win
                plot(ax,tracebp(:, 1),tracebp(:, 2),'r','linew',0.5);
                plot(ax,tracebp(:, 1),tracebp(:, 3),'b','linew',0.5);
                plot(ax,tracebp(:, 1),tracebp(:, 4),'k','linew',0.5);
                
                %%% plot freq. band
                %         plot(ax,[is+1/hihf is+1/lohf],[-0.8*yma -0.8*yma],'k','linew',3)
                if lohf<1 || hihf <1
                    loperd = round(1/lohf);
                    hiperd = round(1/hihf);
                    text(ax,0.02,0.9,strcat(num2str(hiperd),'-',num2str(loperd),{' s'}),'fontsize',...
                        9,'unit',9,'unit','normalized','Horizontalalignment','left');
                else
                    text(ax,0.02,0.9,strcat(num2str(lohf),'-',num2str(hihf),{' Hz'}),'fontsize',9,...
                        'unit','normalized','Horizontalalignment','left');
                end
                if BETALIGN
                    text(ax,0.02,0.75,strcat({'Aligned'}),'fontsize',8,'unit',...
                        'normalized','Horizontalalignment','left');
                    text(ax,0.25,0.1,sprintf('%.2f s',zctime),'fontsize',8,'unit',...
                        'normalized','Horizontalalignment','left');
                else
                    text(ax,0.02,0.75,strcat({'Aligned (poor)'}),'fontsize',8,'unit',...
                        'normalized','Horizontalalignment','left');
                end
                text(ax,0.98,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','right','color','r');
                text(ax,0.98,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','right','color','b');
                text(ax,0.98,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
                    'normalized','Horizontalalignment','right','color','k');
                text(ax,0.35, 0.9, num2str(off13+off13add),'fontsize',8,'color','k','unit','normalized',...
                    'Horizontalalignment','left');
                text(ax,0.25, 0.9, num2str(off12+off12add),'fontsize',8,'color','b','unit','normalized',...
                    'Horizontalalignment','left');
                text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized',...
                    'Horizontalalignment','right');
                
                set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);
                
                if ifreq ~= length(octavef)-1
                    ax.XTickLabel = [];
                else
                    %                 ylabel(ax,'Amplitude','fontsize',9);
                    xlabel(ax,strcat({'Time (s) on '},num2str(date),{', '},DAY,{' '},MO,YEAR),...
                        'fontsize',9);
                end
                %             ax.YTick = [];
                hold(ax,'off');
                
            end
                    
            % %%%%%%%%%%%%%%%%% subplot, amp variation with octave passband %%%%%%%%%%%%%%%%%%%%
            % ax = f.ax(length(octavef)+3);
            % hold(ax,'on');
            % %         ymmat = ymmat/ max(max(ymmat));   % normalized by the global max
            % plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 1)),'r-','linew',1);
            % p1 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 1)),25, 'r','filled');
            % plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 2)),'b-','linew',1);
            % p2 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 2)),25, 'b','filled');
            % plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 3)),'k-','linew',1);
            % p3 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 3)),25, 'k','filled');
            % ymmatmed = median(ymmat, 2);
            % plot(ax,1:size(ymmatmed,1)-1, flipud(ymmatmed(2:end)),'color',[0.6 0.6 0.6],'linew',1);
            % p4 = scatter(ax,1:size(ymmatmed,1), flipud(ymmatmed),25, [0.6 0.6 0.6],'filled');
            % % text(ax,0.02,0.98,sprintf('%s  %s',FLAG,fam),'fontsize',9,'unit','normalized',...
            % %     'Horizontalalignment','left','EdgeColor','k','Margin',0.5);
            % legend(ax,[p1 p2 p3 p4],[stas; 'med '],'location','best');
            % ylabel(ax,'Amplitude','fontsize',9);
            % hold(ax,'off');
            
            %%%%%%%%%%%%%%%%% subplot, power spectral density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first align the broadband trace
            mid = ceil(size(tracebb,1)/2);
            fixlen = 2*sps;
            BETALIGN = 1;
            [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebb(:,2:4)',mid,...
                fixlen,mshiftadd,loffmax,ccmin,iup);
            
            % if a better alignment cannot be achieved, use the LF alignment
            if off12add == mshiftadd+1 && off13add == mshiftadd+1
                off12add = 0;
                off13add = 0;
                BETALIGN = 0;
            end
            % align the records
            istart = sps+1;
            iend = nwinlen+sps;
            tracebblong = [];
            tracebblong(:, 1:2) = tracebb(istart: iend, 1:2);
            tracebblong(:, 3) = tracebb(istart-round(off12add): iend-round(off12add), 3);
            tracebblong(:, 4) = tracebb(istart-round(off13add): iend-round(off13add), 4);
            mid = ceil(size(tracebblong,1)/2);
            
            psdlen = (4:4:16)*sps;
            linesty = {':','-.','--','o-'};
            hlight = 4;
            color=['r','b','k'];
            %     PSDfunc = 'pmtm';
            PSDfunc = 'periodogram';
            Fs = sps;
            
            ax = f.ax(nrow);
            hold(ax, 'on');
            xlim(ax,[1e-1, 1.2*Fs/2]);
            ylim(ax,[-90,0]);
            
            for jj = 1: length(psdlen)
                % cut a segment of different length centered at the main dipole
                x = tracebblong(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
                %     nfft = pow2(nextpow2(size(x,1))-1);
                %     nfft = 1024;
                nfft = size(x,1);
                window = hann(size(x,1));
                nw = 4;
                for ista = 1: nsta
                    if strcmp(PSDfunc, 'pmtm')
                        [psdx,pft] = pmtm(x(:,ista+1),nw,nfft,Fs);
                    elseif strcmp(PSDfunc, 'periodogram')
                        [psdx,pft] = periodogram(x(:,ista+1),window,nfft,Fs);
                    end
                    if jj ~= hlight
                        plot(ax,pft,pow2db(psdx),linesty{jj},'linewidth',0.5,...
                            'color',color(ista),'markers',1.5);
                    else
                        p(ista)=plot(ax,pft,pow2db(psdx),linesty{jj},'linewidth',1,...
                            'color',color(ista),'markers',1.5);
                    end
                end
            end
            % ideally is the reference freq range that the spectrum can be trusted
            plot(ax,[efffreqlo efffreqlo],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
            plot(ax,[efffreqhi efffreqhi],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
            ax.XScale = 'log';
            % text(ax,0.02,0.98,sprintf('%s, %s, %s, Broadband, %d Hz, nfft: %d, Hann',...
            %     FLAG,fam,PSDfunc,sps,nfft),'fontsize',9,'unit','normalized',...
            %     'Horizontalalignment','left');
            xlabel(ax,'Frequency (Hz)');
            ylabel(ax,'PSD (dB/Hz)');
            lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
            hold(ax, 'off');

            %         %save figure
            %         print(f.fig,'-dpdf',strcat(figpath,'/IsoLF_OctaveBP_',fam,'_',FLAG,'_',...
            %               num2str(date),'_',num2str(i),'_',num2str(round(is+idiff/sps)),'s.pdf'));            
        
        end

    end

end





















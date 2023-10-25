%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to identify the LF detections which are featured by an isolated
% arrival. 'isolate' can means the energy ratio between the strongest arrival
% window and its prior window and posterior window. Hopefully we could find
% an appropriate 
%
%   mainly the same as 'identify_merged.m', is for the specific example shown
%   the AGU 2021 abstract!
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/29
% Last modified date:   2021/07/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
clear
close all

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
temppath = strcat(datapath, '/templates/LZBtrio/');
rstpath = strcat(datapath, '/LZBtrio');
figpath = strcat(workpath, '/project2021/TWKBtrio');
hypopath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forproj21');


% combine all fams into an pool, 12 fams in total, final version        
nfampool = ['002';
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
%             '158';      % adding 158 mainly to test
            ];
nfam = size(nfampool, 1);

FLAG = 'TWKB';

% get days, 3-sta trio name, bostock's template name, zero crossing index
% of templates
timoffrot= [
            2003 060;
            2003 061;
            2003 062; % dates of the year, two columns: year Julian_day
            2003 063;
            2003 064;   % LZB only has data from 2003 060-064
            2004 194;
            2004 195;
            2004 196;
            2004 197;
            2004 198;
            2004 199;
            2004 200;
            2004 201;
            2004 202;
            2004 203;
            2005 254;
            2005 255;
            2005 256;
            2005 257;
            2005 258;
            2005 259;
            2005 260;
            2005 261;
            ];
        
nday = size(timoffrot, 1);


%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];

stas=['TWKB'
      'LZB '
      'MGCB'];     % determine the trio and order

%% Important parameters that are consistent with detections
nsta=size(stas,1);         %  number of stations
sps=20;     % samples per second

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensec=16;       % CHECK in detection script for what were used in first-year report
winoffsec=8;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
hi=1.25;
lo=0.5;
npo=2;     % poles, passes of filters
npa=2;

cyclskip = 20;
mshift=14+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmax=4; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnmin=0.35; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?

concentration=1.1; %in seconds; how concentrated is the coherent energy within the window?
cncntr=concentration*sps;   % in samples, 20

%%% load the reference energy ratio and hilbert envelope width from the lfe templates
fname = strcat(temppath, '/temp_engratio.',FLAG,'.sps',num2str(40),'.wlen',num2str(winlensec), ...
            '.f',num2str(lo),'-',num2str(hi));
referat = load(fname);

fname = strcat(temppath, '/temp_hilbwidth.',FLAG,'.sps',num2str(40),'.wlen',num2str(winlensec), ...
            '.f',num2str(lo),'-',num2str(hi));
refhilwid = load(fname);
refhilwid = refhilwid/40*sps;

%%% load detections with or without a distance cutoff
winlensechf = 4;
SUFFIXhf = strcat('up.hf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
                  '.',num2str(xcmaxAVEnmin));
SUFFIXlf = strcat('lf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
                  '.',num2str(xcmaxAVEnmin));

% dcutflag = 0;
dcutflag = 1;

if dcutflag
%     distmaxhf = 10;
    distmaxlf = 12;
%     fnamehf = strcat(hypopath, '/evtloc.pgcfam.pj21.',num2str(distmaxhf),'kmdcutnodou.',SUFFIXhf);
    fnamelf = strcat(hypopath, '/evtloc.lzbfam.pj21.',num2str(distmaxlf),'kmdcutnodou.',SUFFIXlf);
else
%     fnamehf = strcat(hypopath, '/evtloc.pgcfam.pj21.nodcutnodou.',SUFFIXhf);
    fnamelf = strcat(hypopath, '/evtloc.lzbfam.pj21.nodcutnodou.',SUFFIXlf);
end

% allrsthf = load(fnamehf);
allrstlf = load(fnamelf);
% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep + previous 50 cols
%%% UPDATED at 2021/06/23
%%%   n=n+8;
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+12;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)


%% START TO LOOP FOR every day
rstsave = [];

%cycle over each day:
% for nd=9: nday      % num of rows, also num of days
    nd = 9;
    %Which days of data to read?
    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    date = floor(str2double(strcat(YEAR,JDAY)));
    
    rstday = allrstlf(allrstlf(:, 14)==date, :);
    
%     eratprior = rstday(:,21)./rstday(:,22);
%     eratpost = rstday(:,21)./rstday(:,23);
%     eratprior4 = rstday(:,21)./rstday(:,24);
%     eratpost4 = rstday(:,21)./rstday(:,25);
%     eratprior8 = rstday(:,21)./rstday(:,26);
%     eratpost8 = rstday(:,21)./rstday(:,27);
%     eratprior16 = rstday(:,21)./rstday(:,28);
%     eratpost16 = rstday(:,21)./rstday(:,29);
%     % 66 cols, the 43-46th col is the flag of checking by the 4th stations, LZB is least reliable
%     rstday_new = [rstday eratprior eratpost eratprior4 eratpost4 ...
%                   eratprior8 eratpost8 eratprior16 eratpost16];
%     % sort by the energy ratio against prior/posterior 2.5*length_of_strongest_window          
%     [rstday_sort, indori] = sortrows(rstday_new, [59, 60],'descend');
        
    erat(:, 1) = rstday(:,21)./rstday(:,22);
    erat(:, 2) = rstday(:,21)./rstday(:,23);

    % sort by the energy ratio against prior/posterior 2.5*length_of_strongest_window
    % note that the 'erat' has the same index as 'hfbndtime'      
    [erat_sort, indori] = sortrows(erat, [1, 2],'descend');


%     %way 1, choose some fixed thresholds, taken the templates as a reference
%     const = 0.2;
%     indobj = indori(rstday_sort(:,59)>=1.5*rstday_sort(:,60) & ...
%                     rstday_sort(:,60)<=1.5*rstday_sort(:,62) & ...
%                     rstday_sort(:,60)>=const*referat(end, 2) & ...
%                     rstday_sort(:,62)>=const*referat(end, 4) );
%     rstobj = rstday_new(indobj,:);

%     %way 2, just blindly choose the first 5%, but there is a catch, the sorting is in fact
%     %effective to the first sorted col unless there is an same value at the first col. So this first
%     %col first 5% argument mainly is the constaint to the first col only
%     nuse = round(size(rstday_sort,1)* 0.05);
%     indobj = indori(1: nuse);
%     rstobj = rstday_new(indobj,:);
    
    %way 3, get the ?? percentile of col 12 and col 11, based on the actual data, and then use 
    %those as the thresholds 
    perc = 80;
    erat1thr = prctile(erat(:, 1), perc);
    erat2thr = prctile(erat(:, 2), perc);
    indobj = indori(erat_sort(:, 1)>=erat1thr & ...
                    erat_sort(:, 2)>=erat2thr);
    rstobj = rstday(indobj,:);
    
    stasnew=['PGC '
        'SSIB'
        'SILB'
        'KLNB'];  % twkb lzb mgcb
    nstasnew=size(stasnew,1);
    
%     keyboard

    %%% the object might contain detections from different fam, so you need to classify according to
    %%% fam, then read the correct file!
%     for j = 12: size(rstobj,1)
        j = 12;
        famnum = rstobj(j, 13);
        if famnum <= 9
            fam=['00',int2str(famnum)];
        elseif famnum<= 99
            fam=['0',int2str(famnum)];
        else
            fam=int2str(famnum);
        end
        dates = jul2dat(year,jday);
        DAY = int2str(dates(2));
        if dates(1) == 7
            MO = {'Jul. '};
        elseif dates(1) == 3
            MO = {'Mar. '};
        elseif dates(1) == 9
            MO = {'Sep. '};
        end    %%
            
        IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),...
            '.',int2str(npo),int2str(npa),'.ms', int2str(mshift)];
        fname = strcat(rstpath, '/MAPS/pj21traceall_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps','4add');
        
        tracef = load(fname);
        % tracef = [STA1file_new(1:nin_new*winlen,:) STA2file_new(1:nin_new*winlen,2) STA3file_new(1:nin_new*winlen,2) STAopt_new(:, 1:nin_new*winlen)'];
        % time = tracef(:,1);
        % wave1 = tracef(:,2);
        % wave2 = tracef(:,3);
        % wave3 = tracef(:,4);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % maybe consider to load the complete 1-day trace and truncate a segment of the same length and
        % is centered at the main arrival, the advantage is, visually every detection is at the center,
        % and 'findpeaks' would run without any issue when finding the nearest local max and min. But it
        % can also misleads people and myself to think that is the detecting window, which is not true,
        % as the detection might be located in anywhere within the real detecting window (even if j
        % saved the duplicates that are closest to the center of the window
        fname2 = strcat(rstpath, '/MAPS/pj21trace1d_',IDENTIF,'_',num2str(lo),'-',num2str(hi),...
            '_',num2str(winlen/sps),'s',num2str(sps),'sps');
        trace1df = load(fname2);
        % trace1df = [timsSTA' STAopt'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        lofflf = 4;
        ccminlf = 0.35;
        PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
            '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(lofflf),'.ccm', ...
            num2str(ccminlf),'.','80sps');
        fname = strcat(rstpath, '/MAPS/tempfeff_LZB_enc','_',PREFIX);
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
        set(f.ax(8), 'position', [ 0.1, 0.30, 0.85, 0.08]);     % 8-4, for octave filtering
        set(f.ax(7), 'position', [ 0.1, 0.38, 0.85, 0.08]);
        set(f.ax(6), 'position', [ 0.1, 0.46, 0.85, 0.08]);
        set(f.ax(5), 'position', [ 0.1, 0.54, 0.85, 0.08]);
        set(f.ax(4), 'position', [ 0.1, 0.62, 0.85, 0.08]);
        set(f.ax(3), 'position', [ 0.1, 0.70, 0.85, 0.08]);     % for broader band
        set(f.ax(2), 'position', [ 0.1, 0.78, 0.85, 0.08]);     % for LF
        set(f.ax(1), 'position', [ 0.1, 0.86, 0.85, 0.12]);     % for HF
        
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
            [~,ind] = min(abs(rstobj(j, 16)-trace1df(:,1)));
            % here, instead of using 'idiff', we use idiff+ cncntr/2 as a rough estimate to 0-crossing
            ind = ind+cncntr/2;
            trace1d = [];
            trace1d(:,1:2) = trace1df(ind-nwinlen/2: ind+nwinlen/2-1, 1:2);
            % although no pre-alignment, but aligned here according to the results
            off12fec = rstobj(j, 9);
            off13fec = rstobj(j, 10);
            filt12 = filtlf(1);
            filt13 = filtlf(2);
            off12 = off12fec + filt12;
            off13 = off13fec + filt13;
            trace1d(:,3) = trace1df(ind-nwinlen/2 - off12: ind+nwinlen/2-1 - off12, 3);
            trace1d(:,4) = trace1df(ind-nwinlen/2 - off13: ind+nwinlen/2-1 - off13, 4);
            trace = trace1d;
            
            [~,ind] = min(abs(rstobj(j, 16)-trace1dbbf(:,1)));
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
        ocnt = rstobj(j, 15);
        ois = ocnt-winlensec/2;
        oien = ocnt+winlensec/2;
        plot(ax,[ois ois],ax.YLim,'k:','linew',1.5);
        plot(ax,[oien oien],ax.YLim ,'k:','linew',1.5);
        
        %%% plot the seismogram at 3 stations of each win
        plot(ax,trace(:, 1),trace(:, 2),'r','linew',0.5);
        plot(ax,trace(:, 1),trace(:, 3),'b','linew',0.5);
        plot(ax,trace(:, 1),trace(:, 4),'k','linew',0.5);
                
        %%% plot freq. band
%         plot(ax,[is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linew',3)
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
        text(ax,0.75,0.9,sprintf('%.4f',rstobj(j, 17)),'fontsize',8,'unit','normalized',...
             'Horizontalalignment','right');
         % if it is also confirmed to be real by the 4th station
         %         flagconf = [];
         for jj = 1: nstasnew
             if rstobj(j, 43+jj-1) == 1
                 if jj == 1 || jj == 4
                    flagconf(jj) = strcat(stasnew(jj,1));
                 else
                    flagconf(jj) = strcat(stasnew(jj,2));
                 end
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
        ccmin = 0.25;  % 0.4/0.35
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
        %         plot(ax,[is+1/hihf is+1/lohf],[-0.8*yma -0.8*yma],'k','linew',3)
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
%             text(ax,0.75,0.1,sprintf('%.2f/%.2f=%.2f',median(tlag),median(tdura),median(trat)),...
%                 'fontsize',8,'unit','normalized','Horizontalalignment','right');

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
%         octavef = 2.^(3:-1:-5);
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
            % time (j.e., ideally noise or mixture with other arrivals), and the coherency
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
        %         set(f.ax(9), 'position', [ 0.1, 0.22, 0.85, 0.08]);
        set(f.ax(9), 'position', [ 0.1, 0.05, 0.85, 0.21]);     % for spectrum

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
        
%         psdlen = (4:4:16)*sps;
%         hlight = 2;
%         linesty = {':','o-','-.','--'};
        psdlen = 4*sps;
        hlight = 1;
        linesty = {'o-'};
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
                if jj ~= length(hlight)
                    plot(ax,pft,pow2db(psdx),linesty{jj},'linewidth',0.5,...
                        'color',color(ista),'markers',1.5);
                else
                    p(ista)=plot(ax,pft,pow2db(psdx),linesty{jj},'linewidth',1,...
                        'color',color(ista),'markers',1.5);
                end
            end
        end
        % ideally is the reference freq range that the spectrum can be trusted
        if ~isempty(efffreqlo)
            plot(ax,[efffreqlo efffreqlo],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        end
        if ~isempty(efffreqhi)
            plot(ax,[efffreqhi efffreqhi],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        end
        ax.XScale = 'log';
        % text(ax,0.02,0.98,sprintf('%s, %s, %s, Broadband, %d Hz, nfft: %d, Hann',...
        %     FLAG,fam,PSDfunc,sps,nfft),'fontsize',9,'unit','normalized',...
        %     'Horizontalalignment','left');
        xlabel(ax,'Frequency (Hz)');
        ylabel(ax,'PSD (dB/Hz)');
        lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
        hold(ax, 'off');
        
        %save figure
        print(f.fig,'-dpdf',strcat(figpath,'/IsoLF_OctaveBP_PSD_',fam,'_',FLAG,'_',...
            num2str(date),'_',num2str(j),'_',num2str(round(is+idiff/sps)),'s.pdf'));
        
        
        %%%%%%%%%%%%%%%%%%% subplot, plot templates of the same fam %%%%%%%%%%%%%%%%%
%         cla(f.ax(9))
        delete(f.ax(9))
        f.ax(isub) = subplot(nrow,ncol,isub);
            f.ax(isub).FontSize = 8;
            f.ax(isub).Box = 'on';
            grid(f.ax(isub), 'on');
            f.ax(isub).GridLineStyle = '--';
        set(f.ax(9), 'position', [ 0.1, 0.16, 0.85, 0.08]);     % for template

        remake = 0;  % re-make the template or not, 1/0
        dstack = [];
        ccstack = [];
        sps = 40;
        templensec = 60;
        
        disp(fam);
        
        if remake   % if requested templates do not exist, recompute them
            ccmethod = 2; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)
            plflag = 1;
            [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,templensec,ccmethod,plflag);
            ind = [1 2 3];
            %     stack = dstack(ind, :);
            stack = ccstack(ind, :);
            %write into files, NOTE the remade stacks contain all 7 stations
            allstas=['PGC  '
                'SSIB '
                'SILB '
                'LZB  '
                'TWKB '
                'MGCB '
                'KLNB '];
            allnsta = size(allstas,1);
            for ista = 1: allnsta
                fid = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
                    num2str(sps), 'sps_', num2str(templensec), 's_', ...
                    'BBDS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed direct stack, no filter, no norm
                
                fprintf(fid, '%f \n', dstack(ista, :)');
                fclose(fid);
                
                fid = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
                    num2str(sps), 'sps_', num2str(templensec), 's_', ...
                    'BBCCS_', 'opt_Nof_Non_Chao'), 'w+');  % bandpassed cc stack, no filter, no norm
                fprintf(fid, '%f \n', ccstack(ista, :)');
                fclose(fid);
            end
            
        else
            for ista = 1: nsta
                fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                    num2str(templensec), 's_', 'BBDS_', 'opt_Nof_Non_Chao');
                dstack(ista,:) = load(fname);
                fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
                    num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
                ccstack(ista,:) = load(fname);
            end
            stack = ccstack;
        end
        
        templen = size(stack,2);
        mid = templen/2;
        
        wlensec = nwlensec;
        wlen = wlensec*sps;
        tracebb = stack(:, mid-wlen/2-sps: mid+wlen/2+sps-1)';
        
        %%%%%%%% broader band, instead of broadband
        ax = f.ax(9);
        hold(ax,'on');
        
        lohf = 0.5;
        hihf = 6.5;
        npo = 2;
        npa = 2;
        
        TAPER = 'tukey';
        
        if strcmp(TAPER, 'tukey')
            % taper with tukeywin, which is actually a tapered cosine window
            %tapered length is adaptative to frequency, maybe at least longer than one full period length of
            %the lowest frequency
            fractap = round(1/lohf*sps)/size(tracebb,1)*2;
            fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
            w = tukeywin(size(tracebb,1),fractap);
            tracebbtap = w.* tracebb;
        end
        
        tracebpadd = tracebbtap;
        for ista = 1: nsta
            tracebpadd(:,ista) = Bandpass(tracebbtap(:,ista), sps, lohf, hihf, npo, npa, 'butter');
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
        [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd',mid,...
            fixlen,mshiftadd,loffmax,ccmin,iup);
        
        %         % if a better alignment cannot be achieved, use the LF alignment
        %         if off12add == mshiftadd+1 && off13add == mshiftadd+1
        %             off12add = 0;
        %             off13add = 0;
        %             BETALIGN = 0;
        %         end
        
        % align the records
        istart = sps+1;
        iend = wlen+sps;
        tracebp = [];
        tracebp(:, 1) = tracebpadd(istart: iend, 1);
        tracebp(:, 2) = tracebpadd(istart-round(off12add): iend-round(off12add), 2);
        tracebp(:, 3) = tracebpadd(istart-round(off13add): iend-round(off13add), 3);
        
        time = (0:1/sps:wlensec-1/sps)';
        is = round(time(1));  % start time of each win
        ien = round(time(end));   % end time of each win
        
        %in case there is other unpredicted noise, use only the short segment around the dipole
        [maxamppos,indmax] = max(tracebp(max(wlen/2-0.5*fixlen,1): ...
            min(wlen/2+0.5*fixlen-1,wlen),:), [],1);
        [minampneg,indmin] = min(tracebp(max(wlen/2-0.5*fixlen,1): ...
            min(wlen/2+0.5*fixlen-1,wlen),:), [],1);
        if BETALIGN
            %roughly estimate the time of zero-crossing which indicate the origin time of a signal.
            indmax = indmax-1+max(wlen/2-0.5*fixlen,1);  % convert to global index
            indmin = indmin-1+max(wlen/2-0.5*fixlen,1);
            %in case indmin is not smaller than indmax for some detections
            if indmin(1) < indmax(1)
                seg = tracebp(indmin(1): indmax(1), 1);  % for zero-crossing timing, only use the main station
                [~,zc] = min(abs(seg));
                zerocrs1 = find(diff(sign(seg)));
                zc = zc-1+indmin(1);  % convert to global index
            else
                seg = tracebp(indmax(1): indmin(1), 1);  % for zero-crossing timing, only use the main station
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
                    seg = tracebp(indmin(ii)-sps: indmin(ii)-1, ii);
                    tmp = find(diff(sign(seg)),1,'last');
                    if ~isempty(tmp)
                        zclt(ii) = tmp;
                        zclt(ii) = zclt(ii)-1+indmin(ii)-sps;  % convert to global index
                    else
                        zclt(ii) = indmin(ii)-sps;
                    end
                    seg = tracebp(indmax(ii)+1: indmax(ii)+sps, ii);
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
                    seg = tracebp(indmax(ii)-sps: indmax(ii)-1, ii);
                    tmp = find(diff(sign(seg)),1,'last');
                    if ~isempty(tmp)
                        zclt(ii) = tmp;
                        zclt(ii) = zclt(ii)-1+indmax(ii)-sps;  % convert to global index
                    else
                        zclt(ii) = indmax(ii)-sps;
                    end
                    seg = tracebp(indmin(ii)+1: indmin(ii)+sps, ii);
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
        plot(ax,[is+(wlen/2-0.5*fixlen)/sps is+(wlen/2-0.5*fixlen)/sps],ax.YLim,'k--',...
            'linew',1.5);
        plot(ax,[is+(wlen/2+0.5*fixlen-1)/sps is+(wlen/2+0.5*fixlen-1)/sps],ax.YLim,'k--',...
            'linew',1.5);
        
        %%% plot the seismogram at 3 stations of each win
        plot(ax,time,tracebp(:, 1),'r','linew',0.5);
        plot(ax,time,tracebp(:, 2),'b','linew',0.5);
        plot(ax,time,tracebp(:, 3),'k','linew',0.5);
        
        tracetemp = tracebp;
        
        %         text(ax,0.02,0.9,strcat('Broadband'),'fontsize',9,'unit',...
        %              'normalized','Horizontalalignment','left');
        text(ax,0.02,0.9,strcat(num2str(lohf),'-',num2str(hihf),{' Hz'}),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','left');
        if BETALIGN
            text(ax,0.02,0.75,strcat({'Aligned'}),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
%             text(ax,0.25,0.1,sprintf('%.2f s',zctime),'fontsize',8,'unit',...
%                 'normalized','Horizontalalignment','left');
%             %     text(ax,0.75,0.1,sprintf('%.2f/%.2f=%.2f',median(tlag),median(tdura),median(trat)),...
%             %         'fontsize',8,'unit','normalized','Horizontalalignment','right');
%             text(ax,0.75,0.3,sprintf('%.2f/%.2f=%.2f',tlag(1),tdura(1),trat(1)),...
%                 'fontsize',8,'unit','normalized','Horizontalalignment','right','color','r');
%             text(ax,0.75,0.2,sprintf('%.2f/%.2f=%.2f',tlag(2),tdura(2),trat(2)),...
%                 'fontsize',8,'unit','normalized','Horizontalalignment','right','color','b');
%             text(ax,0.75,0.1,sprintf('%.2f/%.2f=%.2f',tlag(3),tdura(3),trat(3)),...
%                 'fontsize',8,'unit','normalized','Horizontalalignment','right','color','k');
            
        else
            text(ax,0.02,0.75,strcat({'Aligned (poor)'}),'fontsize',8,'unit',...
                'normalized','Horizontalalignment','left');
        end
        text(ax,0.02,0.15,sprintf('%s  %s  templates',FLAG,fam),'fontsize',9,'unit','normalized',...
             'Horizontalalignment','left','EdgeColor','k','Margin',0.5);
        text(ax,0.98,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','right','color','r');
        text(ax,0.98,0.75,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','right','color','b');
        text(ax,0.98,0.6,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','right','color','k');
        text(ax,0.35, 0.9, num2str(off13add),'fontsize',8,'color','k','unit','normalized',...
            'Horizontalalignment','left');
        text(ax,0.25, 0.9, num2str(off12add),'fontsize',8,'color','b','unit','normalized',...
            'Horizontalalignment','left');
        text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
            'right');
        
        set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);
        ylabel(ax,'Amplitude','fontsize',9);
        xlabel(ax,'Time (s)','fontsize',9);
        
        % ax.XTickLabel = [];
        %         ax.YTick = [];
        hold(ax,'off');
        
        %save figure
        print(f.fig,'-dpdf',strcat(figpath,'/IsoLF_OctaveBP_temp_',fam,'_',FLAG,'_',...
            num2str(date),'_',num2str(j),'_',num2str(round(is+idiff/sps)),'s.pdf'));


%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     keyboard
%     close all    
% end
        


    
    
    
    
    
    


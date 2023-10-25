% function identify(fam,iup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to identify the LF detections which are featured by an isolated
% arrival. 'isolate' can means the energy ratio between the strongest arrival
% window and its prior window and posterior window. Hopefully we could find
% an appropriate threshold
%
%
% Note 2021/06/10:
%   most oftenly, the signal does not a good SNR and be coherent below 8-16
%   and 16-32 s, so presumably the trace length could be shorter as well,
%   change accordingly
% Note 2021/06/13:
%   And it seems that the prefilter i used for removing the station response
%   is 0.02 0.04 20 40, so the longest untapered period is 25 s, therefore,
%   it doesn't make sense to include pctave bands longer than 25 s, eg, 16-32 s
%
% Note 2021/06/13:
%   Add one more col to the 'rstobj', which is the 34 col of the original result
%   recording the checking flag from the 4th station, KLNB, there are 31-33 cols
%   similarly from PGC trio, but they are too far away to be reliable.
%
% Note 2021/06/17:
%   Is that necessary to use the 40 sps trace, instead of the current 20 sps? Maybe
%   that can give us a better alignment (ie, location) for higher frequency
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/03/15
% Last modified date:   2021/07/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short e   % Set the format to 5-digit floating point
clc
% clear
close all

%% default value for easy debugging
defval('fam', '147');
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
temppath = strcat(datapath, '/templates/LZBtrio/');
rstpath = strcat(datapath, '/LZBtrio');
figpath = strcat(workpath, '/project2021/TWKBtrio');


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


%% START TO LOOP FOR every day
rstsave = [];

for ifam = 1: nfam

    fam = nfampool(ifam, :);
    
    %%% the sequential number of qualified detections from 'way 2' for fam 002
    candidate = goodlfdetections(fam,nday);

%cycle over each day:
for nd=1: nday      % num of rows, also num of days
    
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
    
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),...
             '.',int2str(npo),int2str(npa),'.ms', int2str(mshift)]
    
%%%%%% Uncomment below WHEN using merge_fam_befhypo.m, means removing duplicates BEFORE hypo %%%%%%%        
%     if iup == 1     
%         fname1 = strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '4nsta');   
%         fname2 = strcat(rstpath,'/MAPS/mapallrstNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '3nsta');
%     else
%         fname1 = strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '4nsta');   
%         fname2 = strcat(rstpath,'/MAPS/mapallrstUpNoDou',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
%                         num2str(winlen/sps),'s',num2str(sps),'sps', '3nsta');
%     end
%%%%%% Uncomment above WHEN using merge_fam_befhypo.m, means removing duplicates BEFORE hypo %%%%%%%

%%%%%% Uncomment below WHEN using merge_fam_afthypo.m, means removing duplicates AFTER hypo %%%%%%%        
    if iup == 1     
        fname1 = strcat(rstpath,'/MAPS/pj21mapall',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '4add');   
        fname2 = strcat(rstpath,'/MAPS/pj21mapall',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '3add');
    else
        fname1 = strcat(rstpath,'/MAPS/pj21mapall_up',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '4add');   
        fname2 = strcat(rstpath,'/MAPS/pj21mapall_up',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '3add');
    end
%%%%%% Uncomment above WHEN using merge_fam_befhypo.m, means removing duplicates AFTER hypo %%%%%%%
    
    if isfile(fname1)
        allrst = load(fname1);
        %%% structure of allrst: 30+4*nstanew cols, if 4 new stas, then will be 46 cols
        %%% UPDATED at 2021/03/20
        %%%   1:timswin(n) 2:xmaxSTA12ntmp(n) 3:xmaxSTA13ntmp(n) 4:xcmaxAVEnbang(nin) 5:loopoff(n)
        %%%   6:cumsumtrdiff 7:timswin(n)-winlensec/2+idiff/sps 8:cumsumtrdiff/cumsumtr(winlen)
        %%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
        %%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
        %%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
        %%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
        %%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
        %%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)
        
        %%% additional stations for checking detections
        stasnew=['PGC  '
                 'SSIB '
                 'SILB '
                 'KLNB '];
        if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
            POLSTA(3,:)='KELB ';
            stasnew(4,:)='KELB ';
        else
            POLSTA(3,:)='KLNB ';  % remember to change it back
            stasnew(4,:)='KLNB ';
        end
        nstanew=size(stasnew,1);
    elseif isfile(fname2)
        allrst = load(fname2);
        stasnew=['PGC  '
                 'SSIB '
                 'SILB '];
        nstanew=size(stasnew,1);
    else
        fprintf('Day %s %s of fam %s is skipped because of no detections.\n',YEAR, JDAY, fam);
        continue
    end

    eratprior = allrst(:,9)./allrst(:,10);
    eratpost = allrst(:,9)./allrst(:,11);
    eratprior4 = allrst(:,9)./allrst(:,12);
    eratpost4 = allrst(:,9)./allrst(:,13);
    eratprior8 = allrst(:,9)./allrst(:,14);
    eratpost8 = allrst(:,9)./allrst(:,15);
    eratprior16 = allrst(:,9)./allrst(:,16);
    eratpost16 = allrst(:,9)./allrst(:,17);
    
    date = floor(str2double(strcat(YEAR,JDAY)));
    datemat = date* ones(length(allrst(:,9)),1);
    
    famarr = ceil(str2double(fam))*ones(size(allrst(:,9)));
    
    % 23 cols, the 19-22th col is the flag of checking by the 4th stations, LZB is least reliable
    allrst_new = [datemat allrst(:,1:9) eratprior eratpost eratprior4 eratpost4 ...
                  eratprior8 eratpost8 eratprior16 eratpost16 allrst(:,31:34) famarr];

    % sort by the energy ratio against prior/posterior 2.5*length_of_strongest_window          
    [allrst_sort, indori] = sortrows(allrst_new, [11, 12],'descend');      

%     %way 1, choose some fixed thresholds, taken the templates as a reference
%     const = 0.2;
%     indobj = indori(allrst_sort(:,11)>=1.5*allrst_sort(:,12) & ...
%                     allrst_sort(:,12)<=1.5*allrst_sort(:,14) & ...
%                     allrst_sort(:,12)>=const*referat(end, 2) & ...
%                     allrst_sort(:,14)>=const*referat(end, 4) );
%     rstobj = allrst_new(indobj,:);

%     %way 2, just blindly choose the first 5%, but there is a catch, the sorting is in fact
%     %effective to the first sorted col unless there is an same value at the first col. So this first
%     %col first 5% argument mainly is the constaint to the first col only
%     if nd == 2
%         perc = 0.15;
%     else
%         perc = 0.05;
%     end
%     nuse = round(size(allrst_sort,1)* perc);
%     indobj = indori(1: nuse);
%     rstobj = allrst_new(indobj,:);
    
    %way 3, get the ?? percentile of col 12 and col 11, based on the actual data, and then use 
    %those as the thresholds
    perc = 90;
    col11thr = prctile(allrst_new(:,11), perc);
    col12thr = prctile(allrst_new(:,12), perc);
    indobj = indori(allrst_sort(:,11)>=col11thr & ...
                    allrst_sort(:,12)>=col12thr);
    rstobj = allrst_new(indobj,:);

    
    %%
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
    % as the detection might be located in anywhere within the real detecting window (even if i
    % saved the duplicates that are closest to the center of the window
    fname2 = strcat(rstpath, '/MAPS/pj21trace1d_',IDENTIF,'_',num2str(lo),'-',num2str(hi),...
                       '_',num2str(winlen/sps),'s',num2str(sps),'sps');
    trace1df = load(fname2);
    % trace1df = [timsSTA' STAopt'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load the no filtered trace of each detection window, this had the flaw that it only has the
    % same fixed length as the detection, 16 s, which is not enough for passband 16-32 s
    % Note that all windows are already aligned in LF 
    fname = strcat(rstpath, '/MAPS/pj21traceallbb_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                   num2str(winlen/sps),'s',num2str(sps),'sps','4add');
    tracebbf = load(fname);

    % load the no filtered trace of one day, this would help solve the issue that record length is
    % shorter than that required by the longest period of the filter
    % Note that there is no pre-alignment in this case
    fname = strcat(rstpath, '/MAPS/pj21trace1dbb_',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                   num2str(winlen/sps),'s',num2str(sps),'sps');
    trace1dbbf = load(fname);
            
        
    for i = 1: length(indobj)
        
        %%%%%%%%%%%%%%% PLOT 1, LF detection and its HF counterpart %%%%%%%%%%%%%%%%%%%
        f.fig = figure;
        f.fig.Renderer = 'painters';
        widin = 6;  % maximum width allowed is 8.5 inches
        htin = 6;   % maximum height allowed is 11 inches
        set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
        nrow = 2;
        ncol = 1;
        for isub = 1:nrow*ncol
            f.ax(isub) = subplot(nrow,ncol,isub);
            f.ax(isub).FontSize = 8;
            f.ax(isub).Box = 'on';
            grid(f.ax(isub), 'on');
            f.ax(isub).GridLineStyle = '--';
        end
        
        %%% LF
        ax = f.ax(1);
        hold(ax,'on');
        
        ndet = indobj(i);
        idiff = round((allrst_new(ndet, 8)-allrst_new(ndet, 2)+winlensec/2)*sps);
        trace = tracef(winlen*(ndet-1)+1:winlen*ndet, 1:4);
        
        maxamppos = max(trace(max(idiff-cncntr,1):min(idiff+cncntr,winlen), 2:4), [],1);
        minampneg = min(trace(max(idiff-cncntr,1):min(idiff+cncntr,winlen), 2:4), [],1);    
        ymami = max([maxamppos; -minampneg], [], 1);
        pltscale = 1.2;
        ym = pltscale*max(ymami);

%         %%%%%%%% normalize each trace by its own maximum
%         maxamp = max([maxamppos; -minampneg], [], 1);
%         for ii = 1: 3
%             trace(:, ii+1) = trace(:, ii+1)./maxamp(ii);
%         end
%         yma=2;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        is = trace(1, 1);   % start time of each win
        ien= trace(end, 1);         % end time of each win
        
        axis(ax,[is ien -ym ym]);
        
        %%% patch the main arrival window
        xarea=[is+idiff/sps is+(idiff+cncntr)/sps is+(idiff+cncntr)/sps is+idiff/sps ...
                is+idiff/sps];
        yarea=[-ym -ym ym ym -ym];
        patch(ax,xarea,yarea,'k','Facealpha',0.15,'edgecolor','none');
        
        %%% plot the seismogram at 3 stations of each win
        plot(ax,trace(:, 1),trace(:, 2),'r','linew',1);
        plot(ax,trace(:, 1),trace(:, 3),'b','linew',1);
        plot(ax,trace(:, 1),trace(:, 4),'k','linew',1);
                
        %%% plot freq. band
%         plot(ax,[is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linew',3)
        text(ax,0.02,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'fontsize',9,'unit',...
             'normalized','Horizontalalignment','left');
        text(ax,0.05,0.8,strcat({'Aligned'}),'fontsize',9,'unit',...
             'normalized','Horizontalalignment','left'); 
        text(ax,0.95,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
             'normalized','Horizontalalignment','right','color','r');
        text(ax,0.95,0.8,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
             'normalized','Horizontalalignment','right','color','b');
        text(ax,0.95,0.7,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
             'normalized','Horizontalalignment','right','color','k');
        text(ax,0.02,0.1,sprintf('%s  %s',FLAG,fam),'fontsize',10,'unit','normalized',...
             'Horizontalalignment','left','EdgeColor','k','Margin',0.5);
        off12 = allrst_new(ndet, 3);
        off13 = allrst_new(ndet, 4);
        text(ax,0.3, 0.9, num2str(off13),'fontsize',10,'color','k','unit','normalized',...
             'Horizontalalignment','left');
        text(ax,0.3, 0.8, num2str(off12),'fontsize',10,'color','b','unit','normalized',...
             'Horizontalalignment','left');
        
        set(ax,'XTick',is:4:ien,'fontsize',8);
        ylabel(ax,'Amplitude','fontsize',9);
        %Which days of data to read?
        date = allrst_new(ndet,1);
        year = floor(date/1000);
        jday = floor(date-year*1000);
        dates = jul2dat(year,jday);
        YEAR = int2str(dates(3));
        DAY = int2str(dates(2));
        if dates(1) == 7
            MO = {'Jul. '};
        elseif dates(1) == 3
            MO = {'Mar. '};
        elseif dates(1) == 9
            MO = {'Sep. '};
        end
        xlabel(ax,strcat({'Time (s) on '},num2str(date),{', '},DAY,{' '},MO,YEAR),'fontsize',9);
        %     ax.XMinorTick='on';
        %     ax.YMinorTick='on';
        
        %%%%%%% Below: Find the max of the hilbert transform within the arrival
        %%%%%%%%%%% comparison with 'hilbert' and 'envelope' %%%%%%%%%%%%%%%%%%%
        % using 'abs(hilbert)' should be roughly the same as 'envelope'
        x = trace(:,2:4);
        xmean = mean(x);
        xcentered = bsxfun(@minus,x,xmean);
        % compute envelope amplitude
        xampl = abs(hilbert(xcentered));
        % restore the offset
        wavehilup = bsxfun(@plus,xmean,xampl);
        wavehillo = bsxfun(@minus,xmean,xampl);
        medhilup = median(wavehilup,2); %Median of the 3 hilbert transforms for each sample.
        medhillo = median(wavehillo,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
        [waveenvup, waveenvlo] = envelope(trace(:,2:4));
        
        medenvup = median(waveenvup,2);
        medenvlo = median(waveenvlo,2);
        
        %     plot(ax,time(winlen*(ndet-1)+1:winlen*ndet),medhil,'color',[0.3 0.3 0.3],'linew',0.6); %just to test Hilbert transform
        plot(ax,trace(:, 1),medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
        plot(ax,trace(:, 1),medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);

        [maxhil, loc]=max(medenvup(idiff: min(idiff+cncntr-1, winlen)));
        maxloc=idiff+loc-1;
        icheck=maxloc;
        while icheck<winlen-1 && medenvup(icheck)>0.5*maxhil
            icheck=icheck+1;
        end
        hilwid(i,1)=icheck-maxloc;
        while icheck<winlen-1 && medenvup(icheck)>0.25*maxhil
            icheck=icheck+1;
        end
        hilwid(i,2)=icheck-maxloc;
        while icheck<winlen-1 && medenvup(icheck)>0.125*maxhil
            icheck=icheck+1;
        end
        hilwid(i,3)=icheck-maxloc;
        
        icheck=maxloc;
        while icheck>1 && medenvup(icheck)>0.5*maxhil
            icheck=icheck-1;
        end
        hilwid(i,4)=maxloc-icheck;
        while icheck>1 && medenvup(icheck)>0.25*maxhil
            icheck=icheck-1;
        end
        hilwid(i,5)=maxloc-icheck;
        while icheck>1 && medenvup(icheck)>0.125*maxhil
            icheck=icheck-1;
        end
        hilwid(i,6)=maxloc-icheck;
        
        % This for longer-term averages:
        maxhil=median(medenvup(idiff: min(idiff+cncntr-1, winlen))); %median of arrival window
        maxloc=idiff; %"maxloc" is really start of arrival window
        icheck=maxloc;
        while icheck+cncntr<winlen && median(medenvup(icheck:icheck+cncntr-1))>0.666*maxhil 
            icheck=icheck+1;
        end
        hilwid(i,7)=icheck-maxloc;
        while icheck+cncntr<winlen && median(medenvup(icheck:icheck+cncntr-1))>0.5*maxhil
            icheck=icheck+1;
        end
        hilwid(i,8)=icheck-maxloc;
        while icheck+cncntr<winlen && median(medenvup(icheck:icheck+cncntr-1))>0.25*maxhil
            icheck=icheck+1;
        end
        hilwid(i,9)=icheck-maxloc;
        
        maxloc=idiff+cncntr-1; %"maxloc" is really end of arrival window
        icheck=min(maxloc, winlen);
        while icheck-cncntr>0 && median(medenvup(icheck-cncntr+1:icheck))>0.666*maxhil
            icheck=icheck-1;
        end
        hilwid(i,10)=maxloc-icheck;
        while icheck-cncntr>0 && median(medenvup(icheck-cncntr+1:icheck))>0.5*maxhil 
            icheck=icheck-1;
        end
        hilwid(i,11)=maxloc-icheck;
        while icheck-cncntr>0 && median(medenvup(icheck-cncntr+1:icheck))>0.25*maxhil 
            icheck=icheck-1;
        end
        hilwid(i,12)=maxloc-icheck;
        
        plot(ax, ax.XLim, [0.5*maxhil 0.5*maxhil], '--');
        plot(ax, ax.XLim, [0.25*maxhil 0.25*maxhil], '--');
        plot(ax, ax.XLim, [0.125*maxhil 0.125*maxhil], '--');
        %%%%%%% Above: Find the max of the hilbert transform within the arrival; then walk out from there +/-
        
        %%%%%%% Below: Find the global and nearest local maximum and minimum, compare the amplitude
        %%%%% if use one-day trace
        [~,ind] = min(abs(rstobj(i, 8)-trace1df(:,1)));
        trace1d = trace1df(ind-winlen/2+1: ind+winlen/2, :);
        trace = trace1d;
        idiff = winlen/2;
        %%%%%
          
        maxamppos = max(trace(max(idiff-cncntr,1):min(idiff+cncntr,winlen), 2:4), [],1);
        minampneg = min(trace(max(idiff-cncntr,1):min(idiff+cncntr,winlen), 2:4), [],1);    
        
        for ii = 1: 3
            gmaxloc = find(trace(:, ii+1) == maxamppos(ii));
            [lmaxamp(ii,4:5),loc] = findpeaks(trace(gmaxloc+1: end, ii+1),'NPeaks',2,...
                                                      'MinPeakDistance', 2);
            lmaxloc(ii,4:5) = loc+gmaxloc;                                     
            [lmaxamp(ii,2:-1:1),loc] = findpeaks(flipud(trace(1: gmaxloc-1, ii+1)),'NPeaks',2,...
                                                      'MinPeakDistance', 2);
            lmaxloc(ii,2:-1:1) = gmaxloc-loc;
            lmaxamp(ii,3) = maxamppos(ii);
            lmaxloc(ii,3) = gmaxloc;
            
            gminloc = find(trace(:, ii+1) == minampneg(ii));
            [pk,loc] = findpeaks(-trace(gminloc+1: end, ii+1),'NPeaks',2,...
                                                      'MinPeakDistance', 2);
            lminamp(ii,4:5) = -pk;                                      
            lminloc(ii,4:5) = loc+gminloc;                                     
            [pk,loc] = findpeaks(flipud(-trace(1: gminloc-1, ii+1)),'NPeaks',2,...
                                                      'MinPeakDistance', 2);
            lminamp(ii,2:-1:1) = -pk;
            lminloc(ii,2:-1:1) = gminloc-loc;
            lminamp(ii,3) = minampneg(ii);
            lminloc(ii,3) = gminloc;
        end
        lmaxarat(:, 1:2) = abs(lmaxamp(:, 1:2)./lmaxamp(:,3));
        lmaxarat(:, 3:4) = abs(lmaxamp(:, 4:5)./lmaxamp(:,3));
        lminarat(:, 1:2) = abs(lminamp(:, 1:2)./lminamp(:,3));
        lminarat(:, 3:4) = abs(lminamp(:, 4:5)./lminamp(:,3));
        medlmaxarat = median(lmaxarat, 1);
        medlminarat = median(lminarat, 1);
        lmaxwid(:, 1:2) = lmaxloc(:,3)-lmaxloc(:, 1:2);
        lmaxwid(:, 3:4) = lmaxloc(:,4:5)-lmaxloc(:, 3);
        lminwid(:, 1:2) = lminloc(:,3)-lminloc(:, 1:2);
        lminwid(:, 3:4) = lminloc(:,4:5)-lminloc(:, 3);
        medlmaxwid = median(lmaxwid, 1);
        medlminwid = median(lminwid, 1);
        
        %%%%%%% Above: Find the global and nearest local maximum and minimum, compare the amplitude
        hold(ax,'off');
        
        % the corresponding broadband trace
        tracebb = tracebbf(winlen*(ndet-1)+1:winlen*ndet, 1:4);
        tracebp = tracebb;        
        
        %%%%%% original HF 
        ax = f.ax(2);
        hold(ax,'on');
        
        lohf = 1.25;
        hihf = 6.5;
        npo = 2;
        npa = 2;
        for ista = 1: nsta
            tracebp(:,ista+1) = Bandpass(tracebb(:,ista+1), sps, lohf, hihf, npo, npa, 'butter');
        end
        
%         maxamppos = max(tracehf(max(idiff-cncntr,1):min(idiff+cncntr,winlen), 2:4), [],1);
%         minampneg = min(tracehf(max(idiff-cncntr,1):min(idiff+cncntr,winlen), 2:4), [],1);    
        maxamppos = max(tracebp(:, 2:4), [],1);
        minampneg = min(tracebp(:, 2:4), [],1);
        ymami = max([maxamppos; -minampneg], [], 1);
        ym = pltscale*max(ymami);
        
        %         %%%%%%%% normalize each trace by its own maximum
        %         maxamp = max([maxamppos; -minampneg], [], 1);
        %         for ii = 1: 3
        %             trace(:, ii+1) = trace(:, ii+1)./maxamp(ii);
        %         end
        %         yma=2;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axis(ax,[is ien -ym ym]);
        
        %%% patch the main arrival window
        idiff = round((allrst_new(ndet, 8)-allrst_new(ndet, 2)+winlensec/2)*sps);
        xarea=[is+idiff/sps is+(idiff+cncntr)/sps is+(idiff+cncntr)/sps is+idiff/sps ...
            is+idiff/sps];
        yarea=[-ym -ym ym ym -ym];
        patch(ax,xarea,yarea,'k','Facealpha',0.15,'edgecolor','none');
        
        %%% plot the seismogram at 3 stations of each win
        plot(ax,tracebp(:, 1),tracebp(:, 2),'r','linew',1);
        plot(ax,tracebp(:, 1),tracebp(:, 3),'b','linew',1);
        plot(ax,tracebp(:, 1),tracebp(:, 4),'k','linew',1);
        
        %%% plot freq. band
        %         plot(ax,[is+1/hihf is+1/lohf],[-0.8*yma -0.8*yma],'k','linew',3)
        text(ax,0.02,0.9,strcat(num2str(lohf),'-',num2str(hihf),{' Hz'}),'fontsize',9,'unit',...
            'normalized','Horizontalalignment','left');
        text(ax,0.95,0.9,sprintf('%.4f',ymami(1)),'fontsize',8,'unit',...
             'normalized','Horizontalalignment','right','color','r');
        text(ax,0.95,0.8,sprintf('%.4f',ymami(2)),'fontsize',8,'unit',...
             'normalized','Horizontalalignment','right','color','b');
        text(ax,0.95,0.7,sprintf('%.4f',ymami(3)),'fontsize',8,'unit',...
             'normalized','Horizontalalignment','right','color','k');
        
        %     text(ax,is+3, 0.66*yma, num2str(off13(ndet)),'fontsize',8,'color','k');
        %     text(ax,ien-3, 0.66*yma, num2str(off12(ndet)),'fontsize',8,'color','b');
        
        set(ax,'XTick',is:4:ien,'fontsize',8);
        ylabel(ax,'Amplitude','fontsize',9);
        xlabel(ax,strcat({'Time (s) on '},num2str(date),{', '},DAY,{' '},MO,YEAR),'fontsize',9);
        hold(ax,'off');
%         %%% save figure
%         print(f.fig,'-dpdf',strcat(figpath,'/Isolated_LF_',fam,'_',FLAG,'_',...
%               num2str(date),'_',num2str(i),'_',num2str(round(is+idiff/sps)),'s.pdf'));
        close(f.fig)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        %%%%%%%%%%%%%%%%% PLOT 2, octave-filtered version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f.fig = figure;
        f.fig.Renderer = 'painters';
        widin = 6;  % maximum width allowed is 8.5 inches
        htin = 9.5;   % maximum height allowed is 11 inches
        set(f.fig,'Position',[1*scrsz(3)/15 scrsz(4)/10 widin*res htin*res]);
%         nrow = 11;
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

        if strcmp(WAVE, 'win')
            %%%%% if use trace of the detecting window, pre-aligned in LF
            idiff = round((allrst_new(ndet, 8)-allrst_new(ndet, 2)+winlensec/2)*sps);
            trace = tracef(winlen*(ndet-1)+1:winlen*ndet, 1:4);
            tracebb = tracebbf(winlen*(ndet-1)+1:winlen*ndet, 1:4);
            nwinlen = winlen;
            %%%%%
        elseif strcmp(WAVE, 'day')
            %%%%% if use one-day trace rather than the detecting window
            minfreq = -2;   % lowest oxtave frequency power, ie, lowest freq is 2^-(minfreq)
            nwlensec = max(4*pow2(-minfreq), 32);
            nwinlen = nwlensec*sps;
%             nwinlen = 64*sps;
            [~,ind] = min(abs(rstobj(i, 8)-trace1df(:,1)));
            % here, instead of using 'idiff', we use idiff+ cncntr/2 as a rough estimate to 0-crossing 
            ind = ind+cncntr/2;     
            trace1d = [];
            trace1d(:,1:2) = trace1df(ind-nwinlen/2: ind+nwinlen/2-1, 1:2);
            % although no pre-alignment, but aligned here according to the results
            trace1d(:,3) = trace1df(ind-nwinlen/2 - off12: ind+nwinlen/2-1 - off12, 3);
            trace1d(:,4) = trace1df(ind-nwinlen/2 - off13: ind+nwinlen/2-1 - off13, 4);
            trace = trace1d;
            
            [~,ind] = min(abs(rstobj(i, 8)-trace1dbbf(:,1)));
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

%         %%%%%
%         lo = 0.5;
%         hi = 1.25;
%         npo = 2;
%         npa = 2;
%         for ista = 1: nsta
%             trace(:,ista+1) = Bandpass(tracenof(:,ista+1), sps, lo, hi, npo, npa, 'butter');
%         end
%         %%%%%
        
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
        ym = pltscale*max(ymami);
        axis(ax,[is ien -ym ym]);
        
        %%% patch the main arrival window
        xarea=[is+idiff/sps is+(idiff+cncntr-1)/sps is+(idiff+cncntr-1)/sps is+idiff/sps ...
                is+idiff/sps];
        yarea=[-ym -ym ym ym -ym];
        patch(ax,xarea,yarea,'k','Facealpha',0.15,'edgecolor','none');
        
        % annotate the start and end of the original detecting window 
        temp = tracef(winlen*(ndet-1)+1:winlen*ndet, 1);
        ois = temp(1);
        oien = temp(end);
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
%         text(ax,0.02,0.1,FLAG,'fontsize',10,'unit','normalized','Horizontalalignment','left');
%         text(ax,0.15,0.1,fam,'fontsize',10,'unit','normalized','Horizontalalignment','left');
        text(ax,0.35, 0.9, num2str(off13),'fontsize',8,'color','k','unit','normalized',...
             'Horizontalalignment','left');
        text(ax,0.25, 0.9, num2str(off12),'fontsize',8,'color','b','unit','normalized',...
             'Horizontalalignment','left');
        text(ax,0.75,0.9,sprintf('%.4f',allrst_new(ndet,5)),'fontsize',8,'unit','normalized',...
             'Horizontalalignment','right');
        % if it is also confirmed to be real by the 4th station
%         flagconf = [];     
        for jj = 1: nstanew
            if allrst_new(ndet, 19+jj-1) == 1
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
        plot(ax, ax.XLim, [0.5*maxhil 0.5*maxhil], '--');
        plot(ax, ax.XLim, [0.25*maxhil 0.25*maxhil], '--');
        plot(ax, ax.XLim, [0.125*maxhil 0.125*maxhil], '--');
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
                if jj ~= length(psdlen)
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    close all
end
    keyboard
end
        
%         %%%%%%%%%%%%%%%%% separate PLOT 3, amp variation with octave passband %%%%%%%%%%%%%%%%%%%%
%         f.fig = figure;
%         f.fig.Renderer = 'painters';
%         ax = gca;
%         hold(ax,'on');
%         box(ax,'on');
% %         ymmat = ymmat/ max(max(ymmat));   % normalized by the global max
%         plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 1)),'r-','linew',1);
%         p1 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 1)),25, 'r','filled');
%         plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 2)),'b-','linew',1);
%         p2 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 2)),25, 'b','filled');
%         plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 3)),'k-','linew',1);
%         p3 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 3)),25, 'k','filled');
%         ymmatmed = median(ymmat, 2);
%         plot(ax,1:size(ymmatmed,1)-1, flipud(ymmatmed(2:end)),'color',[0.6 0.6 0.6],'linew',1);
%         p4 = scatter(ax,1:size(ymmatmed,1), flipud(ymmatmed),25, [0.6 0.6 0.6],'filled');
%         text(ax,0.02,0.98,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','left');
%         text(ax,0.2,0.98,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
%         text(ax,0.98,0.98,strcat(num2str(round(is+idiff/sps)),{' s, on '},num2str(date)),...
%              'fontsize',12,'unit','normalized','Horizontalalignment','right');
%         legend(ax,[p1 p2 p3 p4],[stas; 'med '],'location','best');
%         hold(ax,'off');
%         %%% save figure
%         print(f.fig,'-dpdf',strcat(figpath,'/IsoLF_OctaveAmp_',fam,'_',FLAG,'_',...
%               num2str(date),'_',num2str(i),'_',num2str(round(is+idiff/sps)),'s.pdf'));
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         
%         %%%%%%%%%%%%%%%%% separate PLOT 4, power spectral density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %use the broadband trace of the detecting window, so the length is fixed 
%         x = tracebbf(winlen*(ndet-1)+1:winlen*ndet, 1:4);
%         %     PSDfunc = 'pmtm';
%         PSDfunc = 'periodogram';
%         % pmtm param
% %         nfft = pow2(nextpow2(size(x,1))-1);
% %         nfft = 1024;
%         nfft = size(x,1);
%         window = hann(size(x,1));
%         Fs = sps;
%         nw = 4;
%         % filter param
%         npo = 2;
%         npa = 2;
%         lolf = 0.5;
%         lohf = 1.25;
%         hihf = 6.5;
%         
%         f.fig = figure;
%         f.fig.Renderer = 'painters';
%         color=['r','b','k'];
%         ax = gca;
%         hold(ax, 'on');
%         xlim(ax,[1e-2, Fs/2]);
% %         ylim(ax,[-100,0]);
%         for ista = 1: nsta
%             if strcmp(PSDfunc, 'pmtm')
%                 [psdx,pft] = pmtm(x(:,ista+1),nw,nfft,Fs);
%             elseif strcmp(PSDfunc, 'periodogram')
%                 [psdx,pft] = periodogram(x(:,ista+1),window,nfft,Fs);
%             end
%             p(ista)=plot(ax,pft,pow2db(psdx),'o-','linewidth', 1,'color',color(ista),'markers',1.5);
% %             [~,ind] = max(psdx(10:end));
% %             plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(ista));
% %             xhf(:,ista) = Bandpass(tracenof(:,ista), sps, lo, hi, npo, npa, 'butter');
%         end
%         plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         ax.XScale = 'log';
%         hdl = title(sprintf('%s, %s, %s, Broadband, %d Hz, nfft: %d, Hann',...
%                     FLAG,fam,PSDfunc,sps,nfft));
%         movev(hdl,0.02);  
%         text(ax,0.98,0.98,strcat(num2str(round(is+idiff/sps)),{' s, on '},num2str(date)),...
%              'fontsize',12,'unit','normalized','Horizontalalignment','right');
%         xlabel(ax,'Frequency (Hz)');
%         ylabel(ax,'PSD (dB/Hz)');
%         lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
%         grid on
%         box on
%         hold(ax, 'off');
%         %%% save figure
%         print(f.fig,'-dpdf',strcat(figpath,'/IsoLF_BBPSD_',PSDfunc,'_',fam,'_',FLAG,'_',...
%               num2str(date),'_',num2str(i),'_',num2str(round(is+idiff/sps)),'s.pdf'));
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         %%%%%%%%%%%%%%%%% separate PLOT 5, magnitude-squared coherence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % NOT very helpful
%         %coherence, 0 to 1
%         %use the broadband trace of the detecting window, so the length is fixed 
%         x = tracebbf(winlen*(ndet-1)+1:winlen*ndet, 1:4);
%         %     nfft = pow2(nextpow2(size(x,1))-1);
%         nfft = 64;
%         window = hann(nfft);
%         noverlap = nfft*0.5;
%         Fs = sps;
%         f.fig = figure;
%         f.fig.Renderer = 'painters';
%         color=['r','b','k'];
%         ax = gca;
%         hold(ax, 'on');
%         %     ax.XScale = 'log';
%         if strcmp(ax.XScale, 'log')
%             xlim(ax,[1e-2, 1e2]);
%         else
%             xlim(ax,[0, Fs/2]);
%         end
%         ylim(ax,[0,1]);
%         [Cxy(:,1),ft] = mscohere(x(:,1),x(:,2),window,noverlap,nfft,Fs);
%         [Cxy(:,2),~] = mscohere(x(:,1),x(:,3),window,noverlap,nfft,Fs);
%         [Cxy(:,3),~] = mscohere(x(:,2),x(:,3),window,noverlap,nfft,Fs);
%         for ista = 1: nsta
%             p(ista)=plot(ax,ft,Cxy(:,ista),'o-','linewidth', 0.5,'color',color(ista),'markers',1.5);
%         end
%         plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         octavef = 2.^(3:-1:-5);
%         for ifreq = 1: length(octavef)
%             plot(ax,[octavef(ifreq) octavef(ifreq)],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%             if octavef(ifreq) < 1
%                 fstr = strcat(num2str(round(1/octavef(ifreq))), {' s'});
%             else
%                 fstr = strcat(num2str(octavef(ifreq)), {' Hz'});
%             end
%             text(ax,octavef(ifreq),1.02,fstr,'horizontalalignment',...
%                 'center','fontsize',8);
%         end
%         hdl = title(sprintf('%s, %s, Broadband, %d Hz, nfft: %d, noverlap: %d, Hann',...
%                     FLAG,fam,sps,nfft,noverlap));
%         movev(hdl,0.03);  
%         xlabel(ax,'Frequency (Hz)');
%         ylabel(ax,'Magnitude-squared coherence');
%         lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
%         grid on
%         box on
%         hold(ax, 'off');
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         %%%%%%%%%%%%%%%%% separate PLOT 6, phase of cross spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % NOT very helpful
%         %use the broadband trace of the detecting window, so the length is fixed
%         %pre-aligned 
%         x = tracebbf(winlen*(ndet-1)+1:winlen*ndet, 1:4);
%         nfft = 64;
%         window = hann(nfft);
%         noverlap = nfft*0.5;
%         Fs = sps;
%         f.fig = figure;
%         f.fig.Renderer = 'painters';        
%         color=['r','b','k'];
%         ax = gca;
%         hold(ax, 'on');
%         %     ax.XScale = 'log';
%         if strcmp(ax.XScale, 'log')
%             xlim(ax,[1e-2, 1e2]);
%         else
%             xlim(ax,[0, Fs/2]);
%         end
%         ylim(ax,[-1,1]);
%         [Pxy(:,1),ft] = cpsd(x(:,1),x(:,2),window,noverlap,nfft,Fs);
%         [Pxy(:,2),~] = cpsd(x(:,1),x(:,3),window,noverlap,nfft,Fs);
%         [Pxy(:,3),~] = cpsd(x(:,2),x(:,3),window,noverlap,nfft,Fs);
%         for ista = 1: nsta
%             p(ista)=plot(ax,ft,angle(Pxy(:,ista))/pi,'o-','linewidth', 1,'color',...
%                 color(ista),'markers',1.5);
% %         p(ista)=plot(ax,ft,unwrap(angle(Pxy(:,ista)))/pi,'o-','linewidth', 1,'color',...
% %                      color(ista),'markers',1.5);
%         end
%         plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%         hdl = title(sprintf('%s, %s, Broadband, %d Hz, nfft: %d, noverlap: %d, Hann',...
%                     FLAG,fam,sps,nfft,noverlap));
%         movev(hdl,0.02);  
%         xlabel(ax,'Frequency (Hz)');
%         ylabel(ax,'Normalized Phase lag');
%         lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','northeast','fontsize',8);
%         grid on
%         box on
%         hold(ax, 'off');
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
%         % set a tolerance to the width of the hilbert envelope
%         % generally we expect the width of the main lobe of the envelope to be narrow, this is reflected
%         % by col 1-6; we also expect the side lobe amplitude to be small, i.e. long-term average to be
%         % small, this is reflected by col 6-12
%         tolhilwid = 1.5;
%         
%         % set a tolerance to the amplitude ratio between the global max and its two nearst local max
%         % and ratio between the the global min and its two nearst local min (left and right)
%         tolarat = 0.7;
        
%         if sum(hilwid(i,1:6) <= tolhilwid.* refhilwid(end, 1:6)) >= 4 && ...
%            sum(hilwid(i,7:end) <= tolhilwid.* refhilwid(end, 7:end)) >= 4 && ...
%            sum(medlmaxarat <= tolarat) == 4 && ...
%            sum(medlminarat <= tolarat) == 4
%        
%             rstsave = [rstsave; rstobj(i, :)];
%             hold(f.ax(1),'on');
%             plot(f.ax(1),is+(ien-is)*7/8,-0.8*yma, 'kp','MarkerFaceColor','k','MarkerSize',15);
%             hold(f.ax(1),'off');
%         else
%             close(f.fig);
%         end

%     end
        
    
%     %% plot the surrounding LF and HF catalogs of the identified, and checked detections
%     if ~isempty(candidate{nd})
%         figind = setdiff(1:size(rstobj,1), candidate{nd});
%         close(figind);
%       
%         %load the lf catalog of the day
%         lfday = allrst;
%         
%         %load the hf catalog of the day
%         mshifthf=29; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
%         loffhf=2.1; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
%         ccminhf=0.44; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
%         lohi = 1.25;
%         hihf = 6.5;
%         spshf = 40;
%         winlenhf = 4*spshf;
%         IDENTIFhf=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loffhf),'.cc',num2str(ccminhf),...
%              '.',int2str(npo),int2str(npa),'.ms', int2str(mshifthf)];
%         fnamehf = strcat(rstpath,'/MAPS/pj21mapallup',IDENTIFhf,'_',num2str(lohi),'-',...
%                         num2str(hihf),'_',num2str(winlenhf/spshf),'s',num2str(spshf),'sps', '4add');   
%         hfday = load(fnamehf);
%         %%% structure of allrst: 30+4*nstanew cols, if 4 new stas, then will be 46 cols
%         %%% UPDATED at 2021/03/20
%         %%%   1:timswin(n) 2:xmaxSTA12ntmp(n) 3:xmaxSTA13ntmp(n) 4:xcmaxAVEnbang(nin) 5:loopoff(n)
%         %%%   6:cumsumtrdiff 7:timswin(n)-winlensec/2+idiff/sps 8:cumsumtrdiff/cumsumtr(winlen)
%         %%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%         %%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%         %%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%         %%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%         %%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%         %%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)
% 
%         f.fig = figure;
%         f.fig.Renderer = 'painters';
%         widin = 8;  % maximum width allowed is 8.5 inches
%         htin = 8;   % maximum height allowed is 11 inches
%         set(f.fig,'Position',[1*scrsz(3)/15 scrsz(4)/10 widin*res htin*res]);
%         nrow = 2;
%         ncol = 1;
%         for isub = 1:nrow*ncol
%             f.ax(isub) = subplot(nrow,ncol,isub);
%             f.ax(isub).FontSize = 8;
%             f.ax(isub).Box = 'on';
%             grid(f.ax(isub), 'on');
%             f.ax(isub).GridLineStyle = '--';
%         end
%         msize = 15;
%         
%         ax = f.ax(1);
%         hold(ax, 'on');
%         p1=scatter(ax,lfday(:,7),lfday(:,2)*spshf/sps,msize,[0.6 0.6 0.6],'s');  % 12 off hf
%         p2=scatter(ax,lfday(:,7),lfday(:,3)*spshf/sps,msize,[0.6 0.6 0.6],'o','filled');  % 13 off hf
%         scatter(ax,rstobj(candidate{nd},8), rstobj(candidate{nd},3)*spshf/sps,msize*2,...
%                 'bs','filled','markeredgec','k');
%         scatter(ax,rstobj(candidate{nd},8), rstobj(candidate{nd},4)*spshf/sps,msize*2,...
%                 'ro','filled','markeredgec','k');
%         text(ax,0.95,0.1,'LF','horizontalalignment','right','fontsize',12,'unit','normalized');    
%         legend(ax,[p1,p2],{'offset 12','offset 13'},'location','best','fontsize',8);    
%         xlim(ax,[max(min(rstobj(candidate{nd},8))-3600, 0) ...
%                  min(max(rstobj(candidate{nd},8))+3600, 86400)]);
%         ylim(ax,[-40 40])
%         hold(ax, 'off'); 
%         
%         ax = f.ax(2);
%         hold(ax, 'on');
%         p1=scatter(ax,hfday(:,7),hfday(:,2),msize,[0.6 0.6 0.6],'s');  % 12 off hf
%         p2=scatter(ax,hfday(:,7),hfday(:,3),msize,[0.6 0.6 0.6],'o','filled');  % 13 off hf
%         scatter(ax,rstobj(candidate{nd},8), rstobj(candidate{nd},3)*spshf/sps,msize*2,...
%                 'bs','filled','markeredgec','k');
%         scatter(ax,rstobj(candidate{nd},8), rstobj(candidate{nd},4)*spshf/sps,msize*2,...
%                 'ro','filled','markeredgec','k');
%         text(ax,0.95,0.1,'HF','horizontalalignment','right','fontsize',12,'unit','normalized');    
%         xlim(ax,[max(min(rstobj(candidate{nd},8))-3600, 0) ...
%                  min(max(rstobj(candidate{nd},8))+3600, 86400)]);
%         ylim(ax,[-40 40])
%         xlabel(ax,strcat({'Time (s) on '},num2str(date),{', '},DAY,{' '},MO,YEAR),...
%                'fontsize',9);
%         ylabel(ax,'Offset before EFEC (samples)','fontsize',9);  
%         hold(ax, 'off'); 
%     else
%         close all
%     end
%   
%     close all

    
    %%% assemble all checked detections
    rstsave = [rstsave; rstobj(candidate{nd}, :)];

% end
% 
% close all
% end
















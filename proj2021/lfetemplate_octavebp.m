%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to examine the octave bandpass-filtering of the templates
%
% NOTE:
%
% -- Whenever you want to filter a trace, always taper first, this is to avoid
%   the artifact from non-periodicity as filtering in Fourier domain assumes 
%   the incoming signal is periodic which ends at where it starts, and taper
%   would make sure both ends go to 0.  The 'tukeywin' tapered cosine window
%   is a good choice. 
% -- I think maybe 'periodfit' is more useful when the filter bandpass is too
%   narrow so that filtering would be numerically unstable, but you still want
%   to know the contribution from a specific period, then you can directly invert
%   for that period for that small segment of the signal
% -- The phenomenon that the end of the filtered signal goes to 0 but not the
%   start is just the artifact of non-periodicity, although even Frederik does
%   not know why the artifact behaved that way.
% -- Increasing the order of the filter could make the flat portion of the filter
%   bigger, but won't mitigate the artifact, and also will increase the instability
%   of the filter when the passband is too narrow.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/04/28
% Last modified date:   2021/04/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);
% scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
% wid=scrsz(3);       % width
% hite=scrsz(4);       % height

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');

FLAG = 'PGC';
% FLAG = 'TWKB';

disp(FLAG);

if strcmp(FLAG, 'TWKB')
    temppath = strcat(datapath, '/templates/LZBtrio/');
    rstpath = strcat(datapath, '/LZBtrio');
    
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
            ];   
    nfam = size(nfampool,1);

    stas=['TWKB '
          'LZB  '
          'MGCB '];     % determine the trio and order

elseif strcmp(FLAG, 'PGC')    
    temppath = strcat(datapath, '/templates/PGCtrio/');
    rstpath = strcat(datapath, '/PGCtrio');
    
    nfampool = ['002';
                '243';
                '240';
               ];
    nfam = size(nfampool,1);

    stas=['PGC  '
          'SSIB '
          'SILB '];

end
    
% number of used stations
nsta=size(stas,1);         %  number of stations


%% make template or read template
remake = 1;  % re-make the template or not, 1/0
dstack = [];
ccstack = [];
sps = 40;
% templensec = 60;
templensec = 32*4;
% CATA = 'fixed';  % use fixed rots params from Yajun/Allan   
% CATA = 'new';  % use rots computed from new lfe catalog
% CATA = 'old';  % use rots computed from old lfe catalog
% ccbp = [2 8];   % bandpass in CC the raw templates

for ifam = 1: nfam
%     ifam=2;
    fam = nfampool(ifam, :);
    disp(fam);
    
    if remake   % if requested templates do not exist, recompute them
        ccmethod = 2; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)  
        plflag = 0;
        if strcmp(FLAG, 'TWKB')
            ccbp = [2 8];
            [dstack, ccstack] = mk_bbtemp_LZB(fam,sps,templensec,ccmethod,ccbp,plflag);
            ind = [5 4 6];
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

        elseif strcmp(FLAG, 'PGC')  
            if isequal(fam,'002')
                CATA = 'fixed';
                ccbp = [2 8];
                [dstack, ccstack] = mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag);
            else
                CATA = 'new';
                ccbp = [2 8];
                [dstack, ccstack] = mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag);
            end
            ind = [1 2 3];
            stack = ccstack(ind, :);
            
            %write into files, NOTE the remade stacks contain all 7 stations
            allstas=['PGC  '
              'SSIB '
              'SILB '
              'LZB  '
              'TWKB '
              'MGCB '
              'KLNB '];
            nsta = size(allstas,1);

            %write into files
            if isequal(fam,'002') && isequal(CATA,'old')
              suffix = '_catold';
            elseif isequal(fam,'002') && isequal(CATA,'new')
              suffix = '_catnew';
            else
              suffix = [];
            end

            for ista = 1: nsta
              fidds = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBDS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed direct stack, no filter, no norm
              fidccs = fopen(strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBCCS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed cc stack, no filter, no norm
              fprintf(fidds, '%f \n', dstack(ista, :)');
              fclose(fidds);
              fprintf(fidccs, '%f \n', ccstack(ista, :)');
              fclose(fidccs);
            end
        end
        
    else    % if exist, load them directly
        for ista = 1: nsta
          fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
            num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao',suffix);
          ccstack(ista,:) = load(fname);
        end
        stack = ccstack;
    end
    
    templen = size(stack,2);
    mid = templen/2;
    
    wlensec = 120;
    wlen = wlensec*sps;
    tracebb = stack(:, mid-wlen/2-sps: mid+wlen/2+sps-1)';


    %% plot the templates
    %%% plot the 1-step template
    %%% figure 3
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    for ista = 1: nsta
        plot(ax,(1:(wlen+2*sps))/sps,tracebb(:, ista), 'linewidth', 1,'color',...
            color(ista)); hold on
        %     plot(stack(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
%         text(ax,1,0,stas(ista,:),'fontsize',9);
    end
    text(ax,0.95,0.9,'Broadband','fontsize',10,'unit','normalized','Horizontalalignment','right');
    text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
         'right','fontsize',10);
    text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
    text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
    text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Time (s)','fontsize',10);
    ylabel(ax,'Amplitude','fontsize',10);
    xlim(ax,[0 wlensec+2]);
%     ylim(ax,[0 4]);
    ax.Box='on';
    grid(ax, 'on');
    ax.GridLineStyle = '--';
        
    %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
    [envup, envlo] = envelope(tracebb);
    
    medenvup = median(envup,2);
    medenvlo = median(envlo,2);
    
    plot(ax,(1:(wlen+2*sps))/sps,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
    plot(ax,(1:(wlen+2*sps))/sps,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
    hold(ax, 'off');
    
            
    %% filter the templates with octave passbands
    f.fig = figure;
    f.fig.Renderer = 'painters';
    if wlensec <= 64
        widin = 6;
    else
        widin = 8;  % maximum width allowed is 8.5 inches
    end
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 11;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        f.ax(isub).FontSize = 8;
        f.ax(isub).Box = 'on';
        grid(f.ax(isub), 'on');
        f.ax(isub).GridLineStyle = '--';
    end
    
    % reposition
    set(f.ax(11), 'position', [ 0.1, 0.06, 0.85, 0.08]);
    set(f.ax(10), 'position', [ 0.1, 0.14, 0.85, 0.08]);
    set(f.ax(9), 'position', [ 0.1, 0.22, 0.85, 0.08]);
    set(f.ax(8), 'position', [ 0.1, 0.30, 0.85, 0.08]);
    set(f.ax(7), 'position', [ 0.1, 0.38, 0.85, 0.08]);
    set(f.ax(6), 'position', [ 0.1, 0.46, 0.85, 0.08]);
    set(f.ax(5), 'position', [ 0.1, 0.54, 0.85, 0.08]);
    set(f.ax(4), 'position', [ 0.1, 0.62, 0.85, 0.08]);
    set(f.ax(3), 'position', [ 0.1, 0.70, 0.85, 0.08]);
    set(f.ax(2), 'position', [ 0.1, 0.78, 0.85, 0.08]);
    set(f.ax(1), 'position', [ 0.1, 0.86, 0.85, 0.12]);

    %%% LF    
    ax = f.ax(1);
    hold(ax,'on');
    lo = 0.5;
    hi = 1.25;
    npo = 2;
    npa = 2;
    
    % taper with tukeywin, which is actually a tapered cosine window
    %tapered length is adaptative to frequency, maybe at least longer than one full period length of
    %the lowest frequency
    fractap = round(1/lo*sps)/size(tracebb,1)*2;
    fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
    w = tukeywin(size(tracebb,1),fractap);
    tracebbtap = w.* tracebb;

    tracebpadd = tracebbtap;
    for ista = 1: nsta
        tracebpadd(:,ista) = Bandpass(tracebbtap(:,ista), sps, lo, hi, npo, npa, 'butter');
    end

    % re-align the filtered traces, because the alignment was in LF, which may not be true for
    % higher frequency
    mshiftadd = 10;    % maximum allowed shift between 2 traces
    mid = ceil(size(tracebpadd,1)/2);
    fixlen = max(2*sps, 1/lo*sps);
    loffmax = 4;
    ccmin = 0.3;  % 0.4/0.35
    iup = 1;    % times of upsampling
    BETALIGN = 1;
    [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd',mid,...
        fixlen,mshiftadd,loffmax,ccmin,iup);
    
    % if a better alignment cannot be achieved, use the LF alignment
    if off12add == mshiftadd+1 && off13add == mshiftadd+1
        off12add = 0;
        off13add = 0;
        BETALIGN = 0;
    end
    
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
            zc = zc-1+indmin(1);  % convert to global index
        else
            seg = tracebp(indmax(1): indmin(1), 1);  % for zero-crossing timing, only use the main station
            [~,zc] = min(abs(seg));
            zc = zc-1+indmax(1);  % convert to global index
        end
        zctime = (is+zc/sps);
    end
    
    ymami = max([maxamppos; -minampneg], [], 1);
    pltscale = 1.2;
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
    
    [envup, envlo] = envelope(tracebp);
    
    medenvup = median(envup,2);
    medenvlo = median(envlo,2);
    
    plot(ax,time,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
    plot(ax,time,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);

    text(ax,0.02,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'fontsize',9,'unit',...
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
    text(ax,0.02,0.1,sprintf('%s  %s',FLAG,fam),'fontsize',9,'unit','normalized',...
        'Horizontalalignment','left','EdgeColor','k','Margin',0.5);
%     text(ax,0.05,0.1,FLAG,'fontsize',10,'unit','normalized','Horizontalalignment','left');
%     text(ax,0.15,0.1,fam,'fontsize',10,'unit','normalized','Horizontalalignment','left');
    text(ax,0.35, 0.9, num2str(off13add),'fontsize',8,'color','k','unit','normalized',...
        'Horizontalalignment','left');
    text(ax,0.25, 0.9, num2str(off12add),'fontsize',8,'color','b','unit','normalized',...
        'Horizontalalignment','left');
    text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
        'right');
    
    set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);
    ax.XTickLabel = [];
    ylabel(ax,'Amplitude','fontsize',9);
%     ax.YTick = [];
    hold(ax,'off');


    %%%%%% original HF 
    ax = f.ax(2);
    hold(ax,'on');
    lo = 1.25;
    hi = 6.5;
    npo = 2;
    npa = 2;
    
    % taper with tukeywin, which is actually a tapered cosine window
    %tapered length is adaptative to frequency, maybe at least longer than one full period length of
    %the lowest frequency
    fractap = round(1/lo*sps)/size(tracebb,1)*2;
    fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
    w = tukeywin(size(tracebb,1),fractap);
    tracebbtap = w.* tracebb;
    
    tracebpadd = tracebbtap;
    for ista = 1: nsta
        tracebpadd(:,ista) = Bandpass(tracebbtap(:,ista), sps, lo, hi, npo, npa, 'butter');
    end
    
    % re-align the filtered traces, because the alignment was in LF, which may not be true for
    % higher frequency
    %         mshiftadd = 5;    % maximum allowed shift between 2 traces
    mid = ceil(size(tracebpadd,1)/2);
    fixlen = max(2*sps, 1/lo*sps);
    %         loffmax = 4;
    %         ccmin = 0.3;  % 0.4/0.35
    %         iup = 1;    % times of upsampling
    BETALIGN = 1;
    %         [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc(tracebpadd(:,2:4)',mid,fixlen,...
    %                                                                    mshift,loffmax,ccmin);
    [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd',mid,...
        fixlen,mshiftadd,loffmax,ccmin,iup);
    % if a better alignment cannot be achieved, use the LF alignment
    if off12add == mshiftadd+1 && off13add == mshiftadd+1
        off12add = 0;
        off13add = 0;
        BETALIGN = 0;
    end
    
    % align the records
    istart = sps+1;
    iend = wlen+sps;
    tracebp = [];
    tracebp(:, 1) = tracebpadd(istart: iend, 1);
    tracebp(:, 2) = tracebpadd(istart-round(off12add): iend-round(off12add), 2);
    tracebp(:, 3) = tracebpadd(istart-round(off13add): iend-round(off13add), 3);
    
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
            zc = zc-1+indmin(1);  % convert to global index
        else
            seg = tracebp(indmax(1): indmin(1), 1);  % for zero-crossing timing, only use the main station
            [~,zc] = min(abs(seg));
            zc = zc-1+indmax(1);  % convert to global index
        end
        zctime = (is+zc/sps);
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
    
    %%% plot freq. band
    %         plot(ax,[is+1/hihf is+1/lohf],[-0.8*yma -0.8*yma],'k','linew',3)
    text(ax,0.02,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'fontsize',8,'unit',...
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
    text(ax,0.35, 0.9, num2str(off13add),'fontsize',8,'color','k','unit','normalized',...
        'Horizontalalignment','left');
    text(ax,0.25, 0.9, num2str(off12add),'fontsize',8,'color','b','unit','normalized',...
        'Horizontalalignment','left');
    text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
        'right');
    
    set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);
    ax.XTickLabel = [];
    %         ylabel(ax,'Amplitude','fontsize',9);
    %         ax.YTick = [];
    hold(ax,'off');
    
        
    %%%%%%%% broadband, OR, broader-band, 0.5-6.5hz
    ax = f.ax(3);
    hold(ax,'on');
    
    lo = 0.5;
    hi = 6.5;
    npo = 2;
    npa = 2;
    
    % taper with tukeywin, which is actually a tapered cosine window
    %tapered length is adaptative to frequency, maybe at least longer than one full period length of
    %the lowest frequency
    fractap = round(1/lo*sps)/size(tracebb,1)*2;
    fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
    w = tukeywin(size(tracebb,1),fractap);
    tracebbtap = w.* tracebb;

    BBFALG = 0;
    tracebpadd = tracebbtap;
    if ~BBFALG
        for ista = 1: nsta
            tracebpadd(:,ista) = Bandpass(tracebbtap(:,ista), sps, lo, hi, npo, npa, 'butter');
        end
    end
    
    % re-align the filtered traces, because the alignment was in LF, which may not be true for
    % higher frequency
    %         mshift = 10;    % maximum allowed shift between 2 traces
    mid = ceil(size(tracebpadd,1)/2);
    fixlen = max(2*sps, 1/lo*sps);
    %         loffmax = 4;
    %         ccmin = 0.3;  % 0.4/0.35
    %         iup = 1;    % times of upsampling
    BETALIGN = 1;
    [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd',mid,...
        fixlen,mshiftadd,loffmax,ccmin,iup);
    
    % if a better alignment cannot be achieved, use the LF alignment
    if off12add == mshiftadd+1 && off13add == mshiftadd+1
        off12add = 0;
        off13add = 0;
        BETALIGN = 0;
    end
    
    % align the records
    istart = sps+1;
    iend = wlen+sps;
    tracebp = [];
    tracebp(:, 1) = tracebpadd(istart: iend, 1);
    tracebp(:, 2) = tracebpadd(istart-round(off12add): iend-round(off12add), 2);
    tracebp(:, 3) = tracebpadd(istart-round(off13add): iend-round(off13add), 3);
    
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
    
    if BBFALG
        text(ax,0.02,0.9,strcat('Broadband'),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','left');
    else
        text(ax,0.02,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'fontsize',8,'unit',...
            'normalized','Horizontalalignment','left');
    end
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
    text(ax,0.35, 0.9, num2str(off13add),'fontsize',8,'color','k','unit','normalized',...
        'Horizontalalignment','left');
    text(ax,0.25, 0.9, num2str(off12add),'fontsize',8,'color','b','unit','normalized',...
        'Horizontalalignment','left');
    text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
        'right');
    
    set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);
    %         ylabel(ax,'Amplitude','fontsize',9);
    ax.XTickLabel = [];
    %         ax.YTick = [];
    hold(ax,'off');

        
    %%%%%%% HF, octave-passband filtered
%     warning('off');     % turn off warning from 'Bandpass.m'
    minfreq = -5;
    octavef = 2.^(3:-1:minfreq);
    ymmat = zeros(length(octavef),3); %the matrix storing the max/min of the octave-filtered main arrival
    ymmat(1, :) = ymami;    %this is from the broadband above
    for ifreq = 1: length(octavef)-1
        
        ax = f.ax(ifreq+3);
        hold(ax,'on');
        
        lo = octavef(ifreq+1);
        hi = octavef(ifreq);
        npo = 2;    % npo being 2 or 3 is proper
        npa = 2;
        
        % taper with tukeywin, which is actually a tapered cosine window
        %tapered length is adaptative to frequency, maybe at least longer than one full period length of
        %the lowest frequency
        fractap = round(1/lo*sps)/size(tracebb,1)*2;
        fractap = max(fractap,0.1); % if fractap is >=1, n-point von Hann window is returned
        w = tukeywin(size(tracebb,1),fractap);
        tracebbtap = w.* tracebb;

        tracebpadd = tracebbtap;
        if diff(2*[lo hi]/sps)>=0.01  % the bandwidth is moderate
            for ista = 1: nsta
                tracebpadd(:,ista) = Bandpass(tracebbtap(:,ista), sps, lo, hi, npo, npa, 'butter');
            end
        else   % the bandwidth is too narrow, thus the filtering may be less reliable
            nsps = 5;
            [num, denom] = rat(nsps/sps);
            tracebbtapds = zeros(size(tracebbtap,1)*num/denom, nsta);
            tracedsbp = zeros(size(tracebbtap,1)*num/denom, nsta);
            for ista = 1: nsta
                %downsample first
                tracebbtapds(:,ista) = resample(tracebbtap(:,ista),num,denom);   % times = num/denom
                %filter then
                tracedsbp(:,ista) = Bandpass(tracebbtapds(:,ista),nsps,lo,hi,npo,npa,'butter');
                %interpolate finally
                tracebpadd(:,ista) = interp(tracedsbp(:,ista),denom/num);
            end
            text(ax,0.16,0.9,strcat({'narrow passband'}),'fontsize',9,'unit',...
                 'normalized','Horizontalalignment','left');
        end

        % re-align the filtered traces, because the alignment was in LF, which may not be true for
        % higher frequency
        %         mshiftadd = 5;    % maximum allowed shift between 2 traces
        mid = ceil(size(tracebpadd,1)/2);
        fixlen = max(2*sps, 1/lo*sps);
        %         loffmax = 4;
        %         ccmin = 0.3;  % 0.4/0.35
        %         iup = 1;    % times of upsampling
        BETALIGN = 1;
        %         [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc(tracebpadd(:,2:4)',mid,fixlen,...
        %                                                                    mshift,loffmax,ccmin);
        [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebpadd',mid,...
            fixlen,mshiftadd,loffmax,ccmin,iup);
        % if a better alignment cannot be achieved, use the LF alignment
        if off12add == mshiftadd+1 && off13add == mshiftadd+1
            off12add = 0;
            off13add = 0;
            BETALIGN = 0;
        end
        
        % align the records
        istart = sps+1;
        iend = wlen+sps;
        tracebp = [];
        tracebp(:, 1) = tracebpadd(istart: iend, 1);
        tracebp(:, 2) = tracebpadd(istart-round(off12add): iend-round(off12add), 2);
        tracebp(:, 3) = tracebpadd(istart-round(off13add): iend-round(off13add), 3);
        
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
                zc = zc-1+indmin(1);  % convert to global index
            else
                seg = tracebp(indmax(1): indmin(1), 1);  % for zero-crossing timing, only use the main station
                [~,zc] = min(abs(seg));
                zc = zc-1+indmax(1);  % convert to global index
            end
            zctime = (is+zc/sps);
        end
        
        ymami = max([maxamppos; -minampneg], [], 1);
        ym = pltscale*max(ymami);
        ymmat(ifreq+1, :) = ymami;
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
        
        %%% plot freq. band
        %         plot(ax,[is+1/hihf is+1/lohf],[-0.8*yma -0.8*yma],'k','linew',3)
        if lo<1 || hi <1
            lo = round(1/lo);
            hi = round(1/hi);
            text(ax,0.05,0.9,strcat(num2str(hi),'-',num2str(lo),{' s'}),'fontsize',9,'unit',...
                'normalized','Horizontalalignment','left');
        else
            text(ax,0.05,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'fontsize',9,'unit',...
                'normalized','Horizontalalignment','left');
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
        text(ax,0.35, 0.9, num2str(off13add),'fontsize',8,'color','k','unit','normalized',...
            'Horizontalalignment','left');
        text(ax,0.25, 0.9, num2str(off12add),'fontsize',8,'color','b','unit','normalized',...
            'Horizontalalignment','left');
        text(ax,0.75,0.9,sprintf('%.4f',cc),'fontsize',8,'unit','normalized','Horizontalalignment',...
            'right');
        
        set(ax,'XTick',[is is+(ien-is)/4 is+(ien-is)*2/4 is+(ien-is)*3/4 ien],'fontsize',8);

        if ifreq ~= length(octavef)-1
            ax.XTickLabel = [];
        else
            %                 ylabel(ax,'Amplitude','fontsize',9);
            xlabel(ax,strcat({'Time (s)'}),'fontsize',9);
        end
        ax.YTick = [];
        hold(ax,'off');
        
    end
    
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/tempOctaveBP_',fam,'_',FLAG,'_',...
          num2str(sps),'sps_',num2str(templensec), 's','.pdf'));
      
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    ax = gca;
    hold(ax,'on');
    box(ax,'on');
    plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 1)),'r-','linew',1);
    p1 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 1)),25, 'r','filled');
    plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 2)),'b-','linew',1);
    p2 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 2)),25, 'b','filled');
    plot(ax,1:size(ymmat,1)-1, flipud(ymmat(2:end, 3)),'k-','linew',1);
    p3 = scatter(ax,1:size(ymmat,1), flipud(ymmat(:, 3)),25, 'k','filled');
    ymmatmed = median(ymmat, 2);
    plot(ax,1:size(ymmatmed,1)-1, flipud(ymmatmed(2:end)),'color',[0.6 0.6 0.6],'linew',1);
    p4 = scatter(ax,1:size(ymmatmed,1), flipud(ymmatmed),25, [0.6 0.6 0.6],'filled');
    text(ax,0.05,0.9,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','left');
    text(ax,0.2,0.9,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
    legend(ax,[p1 p2 p3 p4],[stas; 'med  '],'location','best');
    hold(ax,'off');
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/tempOctaveAmp_',fam,'_',FLAG,'_',...
          num2str(sps),'sps_',num2str(templensec), 's','.pdf'));
      
    ymmatfam(ifam, :, :) = ymmat;
    

end

%% amplitude variation summary for all fams
if strcmp(FLAG, 'TWKB')
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 6;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 3;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        f.ax(isub).FontSize = 8;
        f.ax(isub).Box = 'on';
        grid(f.ax(isub), 'on');
        f.ax(isub).GridLineStyle = '--';
    end
    
    %     % reposition
    %     set(f.ax(3), 'position', [ 0.1, 0.70, 0.85, 0.08]);
    %     set(f.ax(2), 'position', [ 0.1, 0.78, 0.85, 0.08]);
    %     set(f.ax(1), 'position', [ 0.1, 0.86, 0.85, 0.12]);
    cpool = ['r';'g';'b';'k';'c';'m'];
    sizec = size(cpool,1);
    spool = ['o';'^';'p'];
    
    for ista = 1: nsta
        ax = f.ax(ista);
        hold(ax,'on');
        ampfam(:,:) = ymmatfam(:,:,ista);
        bbamp = ampfam(:,1);
        for ifam = 1: nfam
            if ifam <= sizec
                %plot the bandpassed amplitude
                p(ifam)=plot(ax,1:size(ampfam,2)-1,fliplr(ampfam(ifam,2:end)),cpool(ifam,:),'LineStyle',...
                    '-','linew',0.8,'markers',4,'marker',spool(1,:),'markerfacec',cpool(ifam,:),...
                    'markeredgec',cpool(ifam,:));
                %plot the broadband amplitude
                scatter(ax,9,bbamp(ifam),20,cpool(ifam,:),spool(1,:),'filled');
            elseif ifam <= sizec*2
                p(ifam)=plot(ax,1:size(ampfam,2)-1,fliplr(ampfam(ifam,2:end)),cpool(ifam-sizec,:),...
                    'LineStyle','-','linew',0.8,'markers',4,...
                    'marker',spool(2,:),'markerfacec',cpool(ifam-sizec,:),...
                    'markeredgec',cpool(ifam-sizec,:));
                scatter(ax,9,bbamp(ifam),20,cpool(ifam-sizec,:),spool(2,:),'filled');
            else
                p(ifam)=plot(ax,1:size(ampfam,2)-1,fliplr(ampfam(ifam,2:end)),cpool(ifam-sizec*2,:),...
                    'LineStyle','-','linew',0.8,'markers',6,...
                    'marker',spool(3,:),'markerfacec',cpool(ifam-sizec*2,:),...
                    'markeredgec',cpool(ifam-sizec*2,:));
                scatter(ax,9,bbamp(ifam),25,cpool(ifam-sizec*2,:),spool(3,:),'filled');
            end
            
        end
        text(ax,0.05,0.9,stas(ista,:),'fontsize',12,'unit','normalized','Horizontalalignment','left');
        %     text(ax,0.2,0.9,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
        ylim(ax, [0 0.6]);
        hold(ax,'off');
    end
    legend(f.ax(3),p,nfampool,'location','best');
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/tempOctaveAmp_allfam_',FLAG,'_',...
        num2str(sps),'sps_',num2str(templensec), 's','.pdf'));
    
    %%% median amplitude among 3 stations of all fams
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 6;  % maximum width allowed is 8.5 inches
    htin = 4;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    ax = gca;
    hold(ax,'on');
    ampfammed = [];
    ampfammed(:,:) = median(ymmatfam, 3);
    bbampmed = [];
    bbampmed = ampfammed(:,1);
    for ifam = 1: nfam
        if ifam <= sizec
            %plot the bandpassed amplitude
            p(ifam)=plot(ax,1:size(ampfammed,2)-1,fliplr(ampfammed(ifam,2:end)),cpool(ifam,:),...
                'LineStyle','-','linew',0.8,'markers',4,'marker',spool(1,:),'markerfacec',...
                cpool(ifam,:),'markeredgec',cpool(ifam,:));
            %plot the broadband amplitude
            scatter(ax,9,bbampmed(ifam),20,cpool(ifam,:),spool(1,:),'filled');
        elseif ifam <= sizec*2
            p(ifam)=plot(ax,1:size(ampfammed,2)-1,fliplr(ampfammed(ifam,2:end)),cpool(ifam-sizec,:),...
                'LineStyle','-','linew',0.8,'markers',4,...
                'marker',spool(2,:),'markerfacec',cpool(ifam-sizec,:),...
                'markeredgec',cpool(ifam-sizec,:));
            scatter(ax,9,bbampmed(ifam),20,cpool(ifam-sizec,:),spool(2,:),'filled');
        else
            p(ifam)=plot(ax,1:size(ampfammed,2)-1,fliplr(ampfammed(ifam,2:end)),cpool(ifam-sizec*2,:),...
                'LineStyle','-','linew',0.8,'markers',6,...
                'marker',spool(3,:),'markerfacec',cpool(ifam-sizec*2,:),...
                'markeredgec',cpool(ifam-sizec*2,:));
            scatter(ax,9,bbampmed(ifam),25,cpool(ifam-sizec*2,:),spool(3,:),'filled');
        end
        
    end
    ax.FontSize = 8;
    ax.Box = 'on';
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    ylim(ax, [0 0.4]);
    hold(ax,'off');
    legend(ax,p,nfampool,'location','best');
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/tempOctaveAmp_med_allfam_',FLAG,'_',...
        num2str(sps),'sps_',num2str(templensec), 's','.pdf'));
    
    
    %%% median amplitude of 3 stations for fams at the NW part of region of interest
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 6;  % maximum width allowed is 8.5 inches
    htin = 4;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    ax = gca;
    hold(ax,'on');
    famcheck = [
        '006';
        '017';
        '043';
        '068';
        '099';
        '125';
        '141';
        '144';
        '147';
        ];
    [~,ind] = ismember(famcheck,nfampool,'rows');
    ampfammed = [];
    ampfammed(:,:) = median(ymmatfam(ind,:,:), 3);
    bbampmed = [];
    bbampmed = ampfammed(:,1);
    for i = 1: length(ind)
        if i <= sizec
            %plot the bandpassed amplitude
            p(i)=plot(ax,1:size(ampfammed,2)-1,fliplr(ampfammed(i,2:end)),cpool(i,:),'LineStyle',...
                '-','linew',0.8,'markers',4,'marker',spool(1,:),'markerfacec',cpool(i,:),...
                'markeredgec',cpool(i,:));
            %plot the broadband amplitude
            scatter(ax,9,bbampmed(i),20,cpool(i,:),spool(1,:),'filled');
        elseif i <= sizec*2
            p(i)=plot(ax,1:size(ampfammed,2)-1,fliplr(ampfammed(i,2:end)),cpool(i-sizec,:),...
                'LineStyle','-','linew',0.8,'markers',4,...
                'marker',spool(2,:),'markerfacec',cpool(i-sizec,:),...
                'markeredgec',cpool(i-sizec,:));
            scatter(ax,9,bbampmed(i),20,cpool(i-sizec,:),spool(2,:),'filled');
        else
            p(i)=plot(ax,1:size(ampfammed,2)-1,fliplr(ampfammed(i,2:end)),cpool(i-sizec*2,:),...
                'LineStyle','-','linew',0.8,'markers',6,...
                'marker',spool(3,:),'markerfacec',cpool(i-sizec*2,:),...
                'markeredgec',cpool(i-sizec*2,:));
            scatter(ax,9,bbampmed(i),25,cpool(i-sizec*2,:),spool(3,:),'filled');
        end
        
    end
    ax.FontSize = 8;
    ax.Box = 'on';
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    ylim(ax, [0 0.4]);
    hold(ax,'off');
    legend(ax,p(1: length(ind)),famcheck,'location','best');
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/tempOctaveAmp_med_selfam_',FLAG,'_',...
        num2str(sps),'sps_',num2str(templensec), 's','.pdf'));    
    
end





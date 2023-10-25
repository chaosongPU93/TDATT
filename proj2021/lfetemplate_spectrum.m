% function lfetemplate_spectrum(FLAG,fam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to examine the spectrum of the LFE templates, the frequency response
% of the bandpass filters, and the octave bandpass-filtering of the templates
%
% NOTE:
% -- For the spectrum, we use the 'periodogram', which is the plainest method.
%   we don't use 'pmtm', the multitaper version of PSD, which is the weighted
%   average of the PSDs using different tapers; we also don't use 'pchave',
%   which is a moving-average of PSD using overlapping windows. We think these
%   averaged version of PSD is not suitble for the LFE templates which is 
%   featured by the main dipole in the middle (though this feature might be 
%   due to the stacking of numerous LFEs to increase SNR, which also could 
%   introduce some long-period energy from the interference of different
%   frequencies
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
    
    nfampool = [
                '002';
%                 '243';
%                 '240';
                '047';
               ];
    nfam = size(nfampool,1);

    stas=['PGC  '
          'SSIB '
          'SILB '];

end
    
% number of used stations
nsta=size(stas,1);         %  number of stations


%% make template or read template
remake = 0;  % re-make the template or not, 1/0
dstack = [];
ccstack = [];
sps = 40;
templensec = 60;
% templensec = 32*4;
% CATA = 'fixed';  % use fixed rots params from Yajun/Allan   
CATA = 'new';  % use rots computed from new lfe catalog
% CATA = 'old';  % use rots computed from old lfe catalog
% ccbp = [2 8];   % bandpass in CC the raw templates


for ifam = 1: nfam
%     ifam=2;
    fam = nfampool(ifam, :);
    disp(fam);
    

    if isequal(fam,'002') && isequal(CATA,'old')
      suffix = '_catold';
    elseif isequal(fam,'002') && isequal(CATA,'new')
      suffix = '_catnew';
    else
      suffix = [];
    end
    
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
            allnsta = size(allstas,1);

            %write into files
            for ista = 1: allnsta
              fidds = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
                num2str(sps), 'sps_', num2str(templensec), 's_', ...
                'BBDS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed direct stack, no filter, no norm
              fidccs = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
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
    
    wlensec = 54;
    wlen = wlensec*sps;
    tracebb = stack(:, mid-wlen/2-sps: mid+wlen/2+sps-1)';
    tracebb = detrend(tracebb);

    %% plot the templates
    %%% plot the 1-step template
    %%% figure 3
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    for ista = 1: nsta
        p(ista)=plot(ax,(1:size(tracebb,1))/sps,tracebb(:, ista)+0.5*(ista-1), 'linewidth',... 
            1,'color',color(ista)); hold on
        %     plot(stack(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
        %          color(ista)); hold on
%         text(ax,1,0,stas(ista,:),'fontsize',9);
    end
%     text(ax,0.95,0.9,'Broadband','fontsize',12,'unit','normalized','Horizontalalignment','right');
%     text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'right','fontsize',12);
%     text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
%     text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
%     text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Time (s)','fontsize',12);
    ylabel(ax,'Amplitude','fontsize',12);
    xlim(ax,[0 wlensec]);
%     xlim(ax,[22 38]);
%     ylim(ax,[0 4]);
    ax.Box='on';
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    lgd = legend(ax,p,stas,'location','southeast','fontsize',8);
    titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz'});    %{'nfft: '},num2str(nfft),
    hdl = supertit(ax,titstr,12);
    movev(hdl,0.05);
    
    %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
    [envup, envlo] = envelope(tracebb);
    
    medenvup = median(envup,2);
    medenvlo = median(envlo,2);
    
%     plot(ax,(1:size(tracebb,1))/sps,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
%     plot(ax,(1:size(tracebb,1))/sps,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
    hold(ax, 'off');
    
    %% seeking a sensible measure for the end of coda
    %%%Before doing anything, detrend (including removal of mean), taper, bandpass
    for ista = 1: nsta
      %romve mean, linear trend of template
      tracebbd(:,ista) = detrend(tracebb(:,ista));
      %and taper with tukeywin, which is actually a tapered cosine window
      fractap = 0.15;
      w = tukeywin(size(tracebbd(:,ista),1),fractap);
      tracebbdt(:,ista) = w.* tracebbd(:,ista);
      
      %%% 1.8-4.5 for signal and 1.8-6.0 for template seem like a good pair
      %%% OR, 1.8-5.4 for signal and 1.8-18 for template to keep the most of the template at the high end
      %filter the signal
      hi=18;
      lo=1.8;
      npo=2;
      npa=2;
      trace(:,ista) = Bandpass(tracebbdt(:,ista), sps, lo, hi, npo, npa, 'butter');
    end
    
    %For caution, additional check if the mean is indeed 0, detrend again if not 
    disp(mean(trace))
    trace = detrend(trace);
    
    %%%Calculate the cumulative maximum of absolute amplitude from left & right, and the cumulative 
    %%%standard deviation of amplitude from left & right
    for i = 1: size(trace,1)
      tmp = trace(1:i, :);
      tmp1 = std(tmp,0,1);
      cumstdlt(i,:) = tmp1;
      cummaxlt(i,:) = max(abs(tmp),[],1);
      tmp = trace(i:end, :);
      tmp1 = std(tmp,0,1);
      cumstdrt(i,:) = tmp1;
      cummaxrt(i,:) = max(abs(tmp),[],1);
    end
    figure;
    subplot(411)
    for ista = 1: nsta
      plot((1:size(trace,1))/sps,cumstdlt(:,ista),'color',color(ista)); hold on
    end
    grid on
    title('Cumulative standard deviation of amplitude from left');
    subplot(412)
    for ista = 1: nsta
      plot((1:size(trace,1))/sps,cumstdrt(:,ista),'color',color(ista)); hold on
    end
    grid on
    title('Cumulative standard deviation of amplitude from right');
    subplot(413)
    for ista = 1: nsta
      plot((1:size(trace,1))/sps,cummaxlt(:,ista),'color',color(ista)); hold on
    end
    grid on
    title('Cumulative maximum of absolute amplitude from left');
    subplot(414)
    for ista = 1: nsta
      plot((1:size(trace,1))/sps,cummaxrt(:,ista),'color',color(ista)); hold on
    end
    grid on
    title('Cumulative maximum of absolute amplitude from right'); % this seems most useful in this case
    xlabel('Time (s)');
    
    %%%As found to be most useful from previous plot, do a further differential
    dcummaxrt = diff(cummaxrt);
    figure
    for ista = 1: nsta
      plot((2:size(trace,1))/sps,dcummaxrt(:,ista),'color',color(ista)); hold on
    end
    grid on
    title('Differential cumulative maximum of absolute amplitude from right');
    xlabel('Time (s)');
    ylim([-0.07 0.01]);
    
    %%%Do a running average of the absolute amplitude within a short window. Note that if using the
    %%%raw amplitude, a shorter window may have fewer differences
    reclen = sps/4;   % a quarter sec seems fine
    aa = rectwin(reclen)/reclen;  %normalize to an area of unity 1
    abstrace = abs(trace);  %absolute amplitude
    for ista = 1: nsta
      stave(:,ista) = conv(abstrace(:,ista),aa,'same');
      ltave(ista) = mean([abstrace(1:20*sps,ista); abstrace(40*sps:end,ista)]);
    end
    figure
    plot(trace(:,3),'r-'); hold on
    plot(abstrace(:,3)+0.1,'c-'); hold on
    plot(stave(:,3)+0.2,'k-');
    plot([1 2400],[mad(trace(:,3)) mad(trace(:,3))],'b');
    plot([1 2400],[mad(abstrace(:,3))+0.1 mad(abstrace(:,3))+0.1],'b');
    plot([1 2400],[mad(stave(:,3))+0.2 mad(stave(:,3))+0.2],'b');
    xlim([800 1600]);
    
    figure
    p1=plot(stave(:,1),'r-'); hold on
    p2=plot(stave(:,2),'b-');
    p3=plot(stave(:,3),'k-');
    p4=plot([1 2400],[mad(stave(:,1)) mad(stave(:,1))],'--','color',[.5 .5 .5]);
    plot([1 2400],[mad(stave(:,2)) mad(stave(:,2))],'--','color',[.5 .5 .5]);
    plot([1 2400],[mad(stave(:,3)) mad(stave(:,3))],'--','color',[.5 .5 .5]);
    plot([1438 1438],[0 0.08],'-.','color',[.5 .5 .5]);
    xlim([800 1600]);
    legend([p1,p2,p3,p4],'PGC','SSIB','SILB','MAD');
%     figure
%     plot(stave(:,3)/ltave(:,3));
    
    %% plot the PSD of templates using a variable length around the main dipole
    %first align the broadband trace
    mshiftadd = 5;    % maximum allowed shift between 2 traces
    mid = ceil(size(tracebb,1)/2);
    fixlen = 2*sps;
    loffmax = 4;
    ccmin = 0.3;  % 0.4/0.35
    iup = 1;    % times of upsampling
    BETALIGN = 1;
    [off12add,off13add,cc,iloopoff,loopoff] = constrained_cc_interp(tracebb',mid,...
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
    tracebbali = [];
    tracebbali(:, 1) = tracebb(istart: iend, 1);
    tracebbali(:, 2) = tracebb(istart-round(off12add): iend-round(off12add), 2);
    tracebbali(:, 3) = tracebb(istart-round(off13add): iend-round(off13add), 3);
    mid = ceil(size(tracebbali,1)/2);
    
    %according to the plot, if you able to cover at least 6s after the main dipole, you should be
    %safe to include the effective portion of coda, so that the spectrum reflects the contribution
    %from the main dipole & the coda. It seems that using 16s or 20s centered at the zero-crossing
    %are almost the same for the effective freq range, 0.5-4.5 Hz.
    %So, i think it is fine for the most cases where i show the 16-s sec templates in time & 
    %its spectrum, they are sufficient 
    psdlen = (4:4:20)*sps;
    linesty = {':','-.','--','o-','x-'};
    
    % PSDfunc = 'pmtm';
    PSDfunc = 'periodogram';
    Fs = sps;
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    xlim(ax,[1e-1, 1.2*Fs/2]);
%     xlim(ax,[1e-2, 1e2]);
    ylim(ax,[-90,0]);

    % filter param
    npo = 2;
    npa = 2;
    lo = 1.25;
    hi = 6.5;
    
    for jj = 1: length(psdlen)
        % cut a segment of different length centered at the main dipole
        x = tracebbali(mid-psdlen(jj)/2+1: mid+psdlen(jj)/2, :);
        
        % nfft = pow2(nextpow2(size(x,1))-1);
        % nfft = 1024;
        nfft = size(x,1);
        window = hann(size(x,1));
        nw = 4;
        
        for ista = 1: nsta
            if strcmp(PSDfunc, 'pmtm')
                [psdx,pft] = pmtm(x(:,ista),nw,nfft,Fs);
            elseif strcmp(PSDfunc, 'periodogram')
                [psdx,pft] = periodogram(x(:,ista),window,nfft,Fs);
            end
            if psdlen(jj) ~= 16*sps
                plot(ax,pft,pow2db(psdx),linesty{jj},'linewidth',0.5,...
                    'color',color(ista),'markers',1.5);
            else
                p(ista)=plot(ax,pft,pow2db(psdx),linesty{jj},'linewidth',1,...
                    'color',color(ista),'markers',1.5);
                [~,ind] = max(psdx(10:end));
                plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(ista));
                ftampmaxft(ifam, ista) = pft(ind+9);
                xhf(:,ista) = Bandpass(x(:,ista), sps, lo, hi, npo, npa, 'butter');
            end
        end
    end
    [~,~,~,~,~,HABS2,F,H] = Bandpass(x(:,1), sps, lo, hi, npo, npa, 'butter');
    % when npa = 1, the following 3 ways are the same, but if npa = 2, way 3 is wrong, as the system
    % will filter forward and backward
    %     plot(ax,F,decibel(HABS2),'b-','linew', 2); hold on
    plot(ax,F,pow2db(HABS2),'-','color',[0.6 0.6 0.6],'linew', 2);
    %     plot(ax,F,20*log10(abs(H)),'k--','linew', 1);
    plot(ax,ax.XLim, [-3 -3], '-','color',[0.8 0.8 0.8],'linew', 1.5);
    plot(ax,[lo lo],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    plot(ax,[hi hi],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    ax.XScale = 'log';
    text(ax,0.05,0.3,strcat({'Butterworth: '}),...
        'unit','normalized');
    text(ax,0.05,0.25,strcat(num2str(npo),{' Poles; '},num2str(npa),{' Passes'}),...
        'unit','normalized');
    text(ax,0.05,0.9,PSDfunc,'unit','normalized','fontsize',10);
    %     text(ax,0.95,0.9,'Broadband','unit','normalized','horizontalalignment','right','fontsize',10);
    %     text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
    %          'right','fontsize',10);
    %     text(ax,0.95,0.8,fam,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    %     text(ax,0.95,0.75,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},...
        num2str(min(psdlen)/sps),'-',num2str(max(psdlen)/sps),{' s, '},', Hann');    %{'nfft: '},num2str(nfft),
    hdl = supertit(ax,titstr,12);
    movev(hdl,0.05);
    
    ax.Box='on';
    grid(ax, 'on');
    hold(ax, 'off');
    
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/BBtempPSD_',PSDfunc,'_',fam,'_',FLAG,'_',...
        num2str(sps),'_',num2str(min(psdlen)/sps),'-',num2str(max(psdlen)/sps),'s.pdf'));
      
      
    %% plot the fixed-length PSD of templates & noise before & after

    psdlen = 16*sps;
    linesty = {'o-',':'};
    
    % PSDfunc = 'pmtm';
    PSDfunc = 'periodogram';
    Fs = sps;
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    xlim(ax,[1e-1, 1.2*Fs/2]);
%     xlim(ax,[1e-2, 1e2]);
    ylim(ax,[-90,0]);

    % cut a segment of different length centered at the main dipole
    x = tracebbali(mid-psdlen/2+1: mid+psdlen/2, :);
%     nfft = pow2(nextpow2(size(x,1))-1);
%     nfft = 1024;
    nfft = size(x,1);
    window = hann(size(x,1));
    nw = 4;
    
    for ista = 1: nsta
      if strcmp(PSDfunc, 'pmtm')
        [psdx,pft] = pmtm(x(:,ista),nw,nfft,Fs);
      elseif strcmp(PSDfunc, 'periodogram')
        [psdx,pft] = periodogram(x(:,ista),window,nfft,Fs);
      end
      p(ista)=plot(ax,pft,pow2db(psdx),linesty{1},'linewidth',1,...
        'color',color(ista),'markers',1.5);
      [~,ind] = max(psdx(10:end));
      plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(ista));
    end
    
    % simulate the noise by choose the same length before and after the dipole, average them
    xb = tracebbali(mid-psdlen*3/2+1: mid-psdlen/2, :);
    xa = tracebbali(mid+psdlen/2+1: mid+psdlen*3/2, :);
%     nfft = pow2(nextpow2(size(x,1))-1);
%     nfft = 1024;
    nfft = size(xb,1);
    window = hann(size(xb,1));
    nw = 4;
    
    for ista = 1: nsta
      if strcmp(PSDfunc, 'pmtm')
        [psdxb,pft] = pmtm(xb(:,ista),nw,nfft,Fs);
        [psdxa,~] = pmtm(xa(:,ista),nw,nfft,Fs);
      elseif strcmp(PSDfunc, 'periodogram')
        [psdxb,pft] = periodogram(xb(:,ista),window,nfft,Fs);
        [psdxa,~] = periodogram(xa(:,ista),window,nfft,Fs);
      end
      plot(ax,pft,pow2db((psdxb+psdxa)/2),linesty{2},'linewidth',0.5,...
        'color',color(ista),'markers',1.5);
    end

    ax.XScale = 'log';
    text(ax,0.05,0.9,PSDfunc,'unit','normalized','fontsize',10);
    %     text(ax,0.95,0.9,'Broadband','unit','normalized','horizontalalignment','right','fontsize',10);
    %     text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
    %          'right','fontsize',10);
    %     text(ax,0.95,0.8,fam,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    %     text(ax,0.95,0.75,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},...
        num2str(max(psdlen)/sps),{' s, '},{'nfft: '},num2str(nfft),', Hann');    %
    hdl = supertit(ax,titstr,12);  
%     hdl = title(ax,titstr,'fontsize',12);
    movev(hdl,0.05);
    
    ax.Box='on';
    grid(ax, 'on');
    hold(ax, 'off');
    
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/BBtemp&noiPSD_',PSDfunc,'_',fam,'_',FLAG,'_',...
        num2str(sps),'_',num2str(max(psdlen)/sps),'s.pdf'));    
end

keyboard

    %%% HF-passed PSD
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    xlim(ax,[1e-2, 1e2]);
    ylim(ax,[-100,0]);
    for ista = 1: nsta
        if strcmp(PSDfunc, 'pmtm')
            [psdx,pft] = pmtm(xhf(:,ista),nw,nfft,Fs);
        elseif strcmp(PSDfunc, 'periodogram')
            [psdx,pft] = periodogram(xhf(:,ista),window,nfft,Fs);
        end
        p(ista)=plot(ax,pft,pow2db(psdx),'o-','linewidth', 1,'color',color(ista),'markers',1.5);
        [~,ind] = max(psdx(10:end));
        plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(ista));
        
    end
    plot(ax,[lo lo],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    plot(ax,[hi hi],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    ax.XScale = 'log';
    text(ax,0.05,0.3,strcat({'Butterworth: '}),...
         'unit','normalized');
    text(ax,0.05,0.25,strcat(num2str(npo),{' Poles; '},num2str(npa),{' Passes'}),...
         'unit','normalized');
    text(ax,0.05,0.9,PSDfunc,'unit','normalized','fontsize',10);
    text(ax,0.95,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'unit','normalized',...
         'horizontalalignment','right','fontsize',10);
    text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
         'right','fontsize',10);
    text(ax,0.95,0.8,fam,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    text(ax,0.95,0.75,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    grid on
    box on
    hold(ax, 'off');
    

    % pmtm param
    lo = 0.5;
    hi = 1.25;

    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    xlim(ax,[1e-2, 1e2]);
    ylim(ax,[-100,0]);
    for ista = 1: nsta
        if strcmp(PSDfunc, 'pmtm')
            [psdx,pft] = pmtm(tracebb(:,ista),nw,nfft,Fs);
        elseif strcmp(PSDfunc, 'periodogram')
            [psdx,pft] = periodogram(tracebb(:,ista),window,nfft,Fs);
        end
        p(ista)=plot(ax,pft,pow2db(psdx),'o-','linewidth', 1,'color',color(ista),'markers',1.5);
        [~,ind] = max(psdx(10:end));
        plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(ista));
        
        xhf(:,ista) = Bandpass(tracebb(:,ista), sps, lo, hi, npo, npa, 'butter');
    end
    [~,HABS2,F,H] = Bandpass(tracebb(:,1), sps, lo, hi, npo, npa, 'butter');
    % when npa = 1, the following 3 ways are the same, but if npa = 2, way 3 is wrong, as the system
    % will filter forward and backward
%     plot(ax,F,decibel(HABS2),'b-','linew', 2); hold on
    plot(ax,F,pow2db(HABS2),'-','color',[0.6 0.6 0.6],'linew', 2);
%     plot(ax,F,20*log10(abs(H)),'k--','linew', 1);
    plot(ax,ax.XLim, [-3 -3], '-','color',[0.8 0.8 0.8],'linew', 1.5);
    plot(ax,[lo lo],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    plot(ax,[hi hi],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    ax.XScale = 'log';
    text(ax,0.05,0.3,strcat({'Butterworth: '}),...
         'unit','normalized');
    text(ax,0.05,0.25,strcat(num2str(npo),{' Poles; '},num2str(npa),{' Passes'}),...
         'unit','normalized');
    text(ax,0.05,0.9,PSDfunc,'unit','normalized','fontsize',10);
    text(ax,0.95,0.9,'Broadband','unit','normalized','horizontalalignment','right','fontsize',10);
    text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
         'right','fontsize',10);
    text(ax,0.95,0.8,fam,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    text(ax,0.95,0.75,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    grid on
    box on
    hold(ax, 'off');
    
    %%% LF-passed PSD
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    xlim(ax,[1e-2, 1e2]);
    ylim(ax,[-100,0]);
    for ista = 1: nsta
        if strcmp(PSDfunc, 'pmtm')
            [psdx,pft] = pmtm(xhf(:,ista),nw,nfft,Fs);
        elseif strcmp(PSDfunc, 'periodogram')
            [psdx,pft] = periodogram(xhf(:,ista),window,nfft,Fs);
        end
        p(ista)=plot(ax,pft,pow2db(psdx),'o-','linewidth', 1,'color',color(ista),'markers',1.5);
        [~,ind] = max(psdx(10:end));
        plot(ax,[pft(ind+9) pft(ind+9)], ax.YLim, '--','color',color(ista));
        
    end
    plot(ax,[lo lo],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    plot(ax,[hi hi],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
    ax.XScale = 'log';
    text(ax,0.05,0.3,strcat({'Butterworth: '}),...
         'unit','normalized');
    text(ax,0.05,0.25,strcat(num2str(npo),{' Poles; '},num2str(npa),{' Passes'}),...
         'unit','normalized');
    text(ax,0.05,0.9,PSDfunc,'unit','normalized','fontsize',10);
    text(ax,0.95,0.9,strcat(num2str(lo),'-',num2str(hi),{' Hz'}),'unit','normalized',...
         'horizontalalignment','right','fontsize',10);
    text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
         'right','fontsize',10);
    text(ax,0.95,0.8,fam,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    text(ax,0.95,0.75,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    lgd = legend(ax,p,stas,'location','southwest','fontsize',8);
    grid on
    box on
    hold(ax, 'off');
    
    
    

% end







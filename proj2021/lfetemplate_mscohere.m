%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to examine the magnitude-squared coherence, and cross-spectrum of
% the LFE templates. And do some testing to the phase of cross-spectrum
%
% -- There isn't lots of notes for the 'mscohere', but maybe due to the noise
% in the data, it's better to show the result with linear frequency.
% For the 'cpsd', ideally, if for some frequencies, the time offset between
% traces is the same, then you are supposed to fit a line within these 
% frequencies that would pass through (0,0) via extrapolation. But due to
% noise, the error in the fitted line is large, but it should be at least
% clear to eyes.
%  
% -- When 2*pi*fmin*dt<pi --> dt<1/2fmin --> nshift<Fs/2fmin, the line of phase
% would pass through (0,0), this determines the initial phase shift of the
% starting frequency fmin. when it is larger than that, there will be a 
% periodic effect, so it may not pass through (0,0), depending on how many
% cycles; and meanwhile, the converted lag at the start of frequency would 
% not equal to the true lag, as it has gone through cycles.
% -- Do NOT recommend to use 'unwrap' here
% -- Do NOT recommend to convert the normalized lag to lag in samples with any
% of the method as in 'testcohere.m'
% -- Check 'testcohere.m' for more info
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/04/28
% Last modified date:   2021/05/04
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

% FLAG = 'PGC';
FLAG = 'TWKB';

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
    nfam = length(nfampool);

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
    nfam = length(nfampool);

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
sps = 80;
templensec = 60;
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
    
    wlensec = 60;
    wlen = wlensec*sps;   
    x = stack(:, mid-wlen/2+1: mid+wlen/2)';

%     %% plot the templates
%     %%% plot the 1-step template
%     %%% figure 3
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
%     for ista = 1: nsta
%         plot(ax,(1:wlen)/sps,x(:,ista), 'linewidth', 1,'color',...
%             color(ista)); hold on
%         %     plot(stack(ista, mid-5*sps: mid+5*sps) + 1* ista, 'linewidth', 1,'color',...
%         %          color(ista)); hold on
% %         text(ax,1,0,stas(ista,:),'fontsize',9);
%     end
%     text(ax,0.95,0.9,'Broadband','fontsize',10,'unit','normalized','Horizontalalignment','right');
%     text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'right','fontsize',10);
%     text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
%     text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
%     text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
%     xlabel(ax,'Time (s)','fontsize',10);
%     ylabel(ax,'Amplitude','fontsize',10);
%     xlim(ax,[0 wlensec]);
% %     ylim(ax,[0 4]);
%     ax.Box='on';
%     grid(ax, 'on');
%     ax.GridLineStyle = '--';
%     
%     %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
%     [envup, envlo] = envelope(x);
%     
%     medenvup = median(envup,2);
%     medenvlo = median(envlo,2);
%     
%     plot(ax,(1:wlen)/sps,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
%     plot(ax,(1:wlen)/sps,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
%     hold(ax, 'off');
%     
    
       
    %% obtain magnitude-squared coherence and cross-PSD of templates
    % coherence, 0 to 1
%     nfft = pow2(nextpow2(size(x,1))-1);
    nfft = 2048;
    window = hann(nfft);
    noverlap = nfft*2/4;
    Fs = sps;
    % filter param
    lolf = 0.5;
    lohf = 1.25;
    hihf = 6.5;

    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 8.5;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);    
    nrow = 3;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end
    color=['r','b','k'];
    
    [Cxy(:,1),ft] = mscohere(x(:,1),x(:,2),window,noverlap,nfft,Fs);
    [Cxy(:,2),~] = mscohere(x(:,1),x(:,3),window,noverlap,nfft,Fs);
    [Cxy(:,3),~] = mscohere(x(:,2),x(:,3),window,noverlap,nfft,Fs);
    
    for isub = 1:nrow*ncol
        ax = f.ax(isub);
        hold(ax, 'on');
        %     ax.XScale = 'log';
        if strcmp(ax.XScale, 'log')
            xlim(ax,[1e-2, 1e2]);
        else
            if isub == 1
                xlim(ax,[0, Fs/2]);
            elseif isub == 2
                xlim(ax,[0, 10]);
            else
                xlim(ax,[0, 2]);
            end
            
        end
        ylim(ax,[0,1]);
        for ista = 1: nsta
            p(ista)=plot(ax,ft,Cxy(:,ista),'o-','linewidth', 0.5,'color',color(ista),'markers',1.5);
            %         [~,ind] = max(psdx(10:end));
            %         plot(ax,[ft(ind+9) ft(ind+9)], ax.YLim, '--','color',color(ista));
        end
        %     plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        %     plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        %     plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        octavef = 2.^(3:-1:-5);
        for ifreq = 1: length(octavef)
            plot(ax,[octavef(ifreq) octavef(ifreq)],ax.YLim, '--','color',[0.6 0.6 0.6],'linew',1.5);
            if octavef(ifreq) < 1
                fstr = strcat(num2str(round(1/octavef(ifreq))), {' s'});
            else
                fstr = strcat(num2str(octavef(ifreq)), {' Hz'});
            end
            text(ax,octavef(ifreq),1.02,fstr,'horizontalalignment',...
                'center','fontsize',8);
        end
        %     text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
        %     text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.65,fam,'fontsize',10,'unit','normalized','Horizontalalignment','left');
        %     text(ax,0.05,0.6,FLAG,'fontsize',10,'unit','normalized','Horizontalalignment','left');
        xlabel(ax,'Frequency (Hz)');
        ylabel(ax,'Magnitude-squared coherence');
        lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
        grid(ax, 'on');
        ax.Box = 'on';
        hold(ax, 'off');
    end
    titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},num2str(wlensec),...
                    {' s, '},{'nfft: '},num2str(nfft),{', noverlap: '},num2str(noverlap),', Hann');
    hdl = title(f.ax(1),titstr,'fontsize',12);
    movev(hdl,0.05);

    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/BBtempmscohe_',fam,'_',FLAG,'_',...
          num2str(sps),'_',num2str(wlensec),'s.pdf'));  


    
    %% phase of cross spectrum in normalized phase lag
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 8.5;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 3;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
    end    
    color=['r','b','k'];
    
    [Pxy(:,1),ft] = cpsd(x(:,1),x(:,2),window,noverlap,nfft,Fs);
    [Pxy(:,2),~] = cpsd(x(:,1),x(:,3),window,noverlap,nfft,Fs);
    [Pxy(:,3),~] = cpsd(x(:,2),x(:,3),window,noverlap,nfft,Fs);
    
    for isub = 1:nrow*ncol
        ax = f.ax(isub);
        hold(ax, 'on');
        %     ax.XScale = 'log';
        if strcmp(ax.XScale, 'log')
            xlim(ax,[1e-2, 1e2]);
        else
            if isub == 1
                xlim(ax,[0, Fs/2]);
            elseif isub == 2
                xlim(ax,[0, 10]);
            else
                xlim(ax,[0, 2]);
            end
            
        end
        ylim(ax,[-1,1]);
        for ista = 1: nsta
            p(ista)=plot(ax,ft,angle(Pxy(:,ista))/pi,'o-','linewidth', 1,'color',...
                color(ista),'markers',1.5);
            %         p(ista)=plot(ax,ft,unwrap(angle(Pxy(:,ista)))/pi,'o-','linewidth', 1,'color',...
            %                      color(ista),'markers',1.5);
        end
        plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        %     text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
        %     text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        %     text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
        %          'left','fontsize',10);
        xlabel(ax,'Frequency (Hz)');
        ylabel(ax,'Normalized Phase lag');
        lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
        grid(ax, 'on');
        ax.Box = 'on';
        hold(ax, 'off');
        
    end
    titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},num2str(wlensec),...
                    {' s, '},{'nfft: '},num2str(nfft),{', noverlap: '},num2str(noverlap),', Hann');
    hdl = title(f.ax(1),titstr,'fontsize',12);
    movev(hdl,0.05);

    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/BBtempCPSD_',fam,'_',FLAG,'_',...
          num2str(sps),'_',num2str(wlensec),'s.pdf'));    

      
    %% linear fitting to the phase of cross-spectrum
    famcheck = ['043';
                '068';
                '147';
                '141';
                '099';
                '006';
                '125';
                '017';
                '144';];
    fi = 1;

    if ismember(fam,famcheck,'rows')
        switch fam
            case '043'
                ran1 = [0.5 1.25];
                ran2 = [1.25 6];
            case '141'
                ran1 = [0.5 1.25];
                ran2 = [1.25 5];
            case '144'
                ran1 = [0.3 1.25];
                ran2 = [1.25 6];
            case '099'
                ran1 = [0.5 1.25];
                ran2 = [2.5 5];
            case '068'
                ran1 = [0.5 1];
                ran2 = [1 5];
            case '125'
                ran1 = [0.5 1];
                ran2 = [3 5];
            case '147'
                ran1 = [0.5 1.25];
                ran2 = [1.25 3.5];
            case '017'
                ran1 = [0.3 1.25];
                ran2 = [1.25 4];
            case '006'
                ran1 = [0.5 1.25];
                ran2 = [2.5 5.5];
        end
        [~, indst1] = min(abs(ft-ran1(1)));
        [~, inded1] = min(abs(ft-ran1(end)));
        [~, indst2] = min(abs(ft-ran2(1)));
        [~, inded2] = min(abs(ft-ran2(end)));
        
        %create fit object with free constraints
        fttpfree = fittype( @(a,b,x) a*x+b);
        %create fit object with fix constraints
        b = 0;
        fttpfix = fittype( @(a,x) a*x+b);
        
        f.fig = figure;
        f.fig.Renderer = 'painters';
        widin = 8.5;  % maximum width allowed is 8.5 inches
        htin = 9;   % maximum height allowed is 11 inches
        set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
        nrow = nsta;
        ncol = 1;
        for isub = 1:nrow*ncol
            f.ax(isub) = subplot(nrow,ncol,isub);
        end

        lgdstr = {'TWKB-LZB','TWKB-MGCB','LZB-MGCB'};
        for ista = 1: nsta
            %first fit the slope and intercept of the higher-freq range
            x2 = ft(indst2: inded2);
            %     phamed = median(angle(Pxy)/pi,2);
            %     y2 = phamed(indst2: inded2);
            y2 = angle(Pxy(indst2: inded2, ista))/pi;
            
            %linear robust least square
            [fitobj,gof,~] = fit(x2,y2,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
            %output fit parameters
            coef = coeffvalues(fitobj);
            slope2 = coef(1);
            intcpt2 = coef(2);
            
            ax = f.ax(ista);
            hold(ax, 'on');
            grid(ax, 'on');
            ax.Box = 'on';
            xlim(ax,[0, 8]);
            ylim(ax,[-1,1]);
            plot(ax,[ran1(1) ran1(1)],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
            plot(ax,[ran1(end) ran1(end)],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
            plot(ax,[ran2(1) ran2(1)],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
            plot(ax,[ran2(end) ran2(end)],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
            p1 = plot(ax,ft,angle(Pxy(:,ista))/pi,'o-','color',[0.3 0.3 0.3],'linew',1,'markers',1.5);
            yfit2 = linefcn(x2,slope2,intcpt2);
            p2 = plot(ax,x2,yfit2,'-','color','r','linew',1.5);
            
            %then fit the slope of the lower-freq range
            x1 = ft(indst1: inded1);
            %     y1 = phamed(indst1: inded1);
            y1 = angle(Pxy(indst1: inded1, ista))/pi;
            
            %linear robust least square
            [fitobj,gof,~] = fit(x1,y1,fttpfix,'Robust','Bisquare','StartPoint',1);
            %output fit parameters
            coef = coeffvalues(fitobj);
            slope1 = coef(1);
            
            yfit1 = linefcn(x1,slope1,b);
            p3 = plot(ax,x1,yfit1,'-','color','b','linew',1.5);
            xlabel(ax,'Frequency (Hz)');
            ylabel(ax,'Normalized Phase lag');
            lgd = legend(ax,p1,lgdstr{ista},'location','north','fontsize',8);
            hold(ax, 'off');
            
            %%%the theoretical prediction of shifted time or samples based on the fitting, check the NOTE
            %%%'phase of the cross-spectrum' for detailed derivation
            %%%HERE, the sign of tshift or nshift is positive if the trace 2 in cpsd(tr1,tr2) is
            %%%lagging, negative if trace 2 is leading
            %%%compared the sign convention with the filtering effect, e.g., xmaxstack12ntmplf, if it is
            %%%positive, then it means the 2nd station is leading, otherwise the 1st station is leading.
            %%%Take fam 144 as an example, for fi=1, the tshift and nshift are positive for all station
            %%%pairs, then it means, TWKB is leading LZB, TWKB is leading MGCB, and LZB is leading MGCB.
            %%%Meanwhile, the filtering effect shows that in LF, xmaxstack12ntmplf<0, xmaxstack13ntmplf<0,
            %%%xmaxstack32ntmplf>0, i.e., the derived relationships, leading or lagging between stations
            %%%are consistent.
            
            %give a value of frequency
            if fi >= ran1(1) && fi <= ran1(end)
                tshift(ifam, ista) = slope1/2;
                nshift(ifam,ista) = tshift(ifam,ista) *Fs;
                
            elseif fi > ran2(1) && fi <= ran2(end)
                tshift(ifam,ista) = slope2/2+intcpt2/2/fi;
                nshift(ifam,ista) = tshift(ifam,ista) *Fs;
                
            end
            
        end
        titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},{'nfft: '},...
                        num2str(nfft),{', noverlap: '},num2str(noverlap),', Hann');
        hdl = title(f.ax(1),titstr,'fontsize',12);
        movev(hdl,0.05);

        %by comparing with the sign convention of the CC in detection and filtering effect, we need to
        %modify the sign here to be consistent.
        tshift(ifam,:) = -tshift(ifam,:);
        nshift(ifam,:) = -nshift(ifam,:);
        
        for ista = 1: nsta
            disp(nshift(ifam,ista));
            if nshift(ifam,ista) <=0
                fprintf('2nd station in pair %s is lagging by %d samples \n',lgdstr{ista},...
                    abs(round(nshift(ifam,ista))));
            else
                fprintf('1st station in pair %s is lagging by %d samples\n', lgdstr{ista},...
                    abs(round(nshift(ifam,ista))));
            end
        end
        
    else
        tshift(ifam,:) = nan(1,nsta);
        nshift(ifam,:) = nan(1,nsta);
    end

    
end
    %%
%     %%% phase of cross spectrum in lag samples
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
% %     ax.XScale = 'log';
%     if strcmp(ax.XScale, 'log')
%         xlim(ax,[1e-2, 1e2]);
%     else
%         xlim(ax,[0, 10]);
%     end
%     ylim(ax,[-10,10]);
%     [Pxy(:,1),ft] = cpsd(x(:,1),x(:,2),window,noverlap,nfft,Fs);
%     [Pxy(:,2),~] = cpsd(x(:,1),x(:,3),window,noverlap,nfft,Fs);
%     [Pxy(:,3),~] = cpsd(x(:,2),x(:,3),window,noverlap,nfft,Fs);
%     for ista = 1: nsta
% %         p(ista)=plot(ax,ft,angle(Pxy(:,ista))/pi/2*Fs./ft,'o-','linewidth', 1,'color',...
% %                      color(ista),'markers',1.5);
%         lag = diff(unwrap(angle(Pxy(:,ista)))/pi)/2*Fs./(ft(2)-ft(1));
%         p(ista)=plot(ax,ft,[0; lag],'o-','linewidth', 1,'color',...
%                     color(ista),'markers',1.5);
%     end
%     plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);    
%     plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
%     text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10); 
%     xlabel(ax,'Frequency (Hz)');
%     ylabel(ax,'Lag (sample)');
%     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
%     grid on
%     box on
%     hold(ax, 'off');
    





%     %% TEST, make a copy of templates but with some shift
%     nshift = 10;
%     
%     xori = x(1:end-nshift,1);
%     xsft = x(1+nshift:end,1);
%     t = (0:wlen-1-nshift)/Fs;
% 
%     %%% plot two traces
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
%         plot(ax,t,xori, 'linewidth', 1,'color',...
%             color(1)); hold on
%         plot(ax,t,xsft, 'linewidth', 1,'color',...
%             color(2));
%     text(ax,0.95,0.9,'Broadband','fontsize',10,'unit','normalized','Horizontalalignment','right');
%     text(ax,0.95,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'right','fontsize',10);
%     text(ax,0.03,0.92,'a','FontSize',9,'unit','normalized','EdgeColor','k','Margin',2);
%     text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
%     text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
%     xlabel(ax,'Time (s)','fontsize',10);
%     ylabel(ax,'Amplitude','fontsize',10);
%     xlim(ax,[0 wlensec]);
% %     ylim(ax,[0 4]);
%     ax.Box='on';
%     grid(ax, 'on');
%     ax.GridLineStyle = '--';
%     
%     
%     %%% obtain magnitude-squared coherence and cross-PSD of templates
%     % coherence, 0 to 1
% %     nfft = pow2(nextpow2(size(x,1))-1);
%     nfft = 512;
%     window = hann(nfft);
%     noverlap = nfft*2/4;
%     Fs = sps;
%     % filter param
%     lolf = 0.5;
%     lohf = 1.25;
%     hihf = 6.5;
% 
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
% %     xlim(ax,[1e-2, 1e2]);
%     xlim(ax,[0, Fs/2]);
%     ylim(ax,[0,1.2]);
%     [Cxy,ft] = mscohere(xori,xsft,window,noverlap,nfft,Fs);
%         p(ista)=plot(ax,ft,Cxy,'o-','linewidth', 1,'color',color(ista),'markers',1.5);
% %         [~,ind] = max(psdx(10:end));
% %         plot(ax,[ft(ind+9) ft(ind+9)], ax.YLim, '--','color',color(ista));
%     plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);    
%     plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
% %     ax.XScale = 'log';
%     text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
%     text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     xlabel(ax,'Frequency (Hz)');
%     ylabel(ax,'Magnitude-squared coherence');
% %     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
%     grid on
%     box on
%     hold(ax, 'off');
% 
%     
%     %%% phase of cross spectrum in normalized phase lag
%     %%% the slope of the normalized phase lag = nshift / F_Nyquist, it is most reliable, even though
%     %%% sometimes the best-fit line does not pass through 0 hz and 0 lag
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
% %     xlim(ax,[1e-2, 1e2]);
%     xlim(ax,[0, Fs/2]);
% %     ylim(ax,[-1.2,1.2]);
%     [Pxy,ft] = cpsd(xori,xsft,window,noverlap,nfft,Fs);
%         p(ista)=plot(ax,ft,unwrap(angle(Pxy))/pi,'o-','linewidth', 1,'color',...
%                      color(ista),'markers',1.5);
% %         p(ista)=plot(ax,ft,angle(Pxy)/pi,'o-','linewidth', 1,'color',...
% %                      color(ista),'markers',1.5);
%     plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);    
%     plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
% %     ax.XScale = 'log';
%     text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
%     text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     xlabel(ax,'Frequency (Hz)');
%     ylabel(ax,'Normalized Phase lag');
% %     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
%     grid on
%     box on
%     hold(ax, 'off');
%     
%     %%% phase of cross spectrum in lag samples
%     %%% the lag samples would converge to the correct value, but it is not reliable all the time, for
%     %%% example, if nshift = 15 samples, Fs=80, it would fail
%     f.fig = figure;
%     f.fig.Renderer = 'painters';
%     color=['r','b','k'];
%     ax = gca;
%     hold(ax, 'on');
% %     xlim(ax,[1e-2, 1e2]);
%     xlim(ax,[0, Fs/2]);
% %     ylim(ax,[-1.5*nshift,1.5*nshift]);
%     [Pxy,ft] = cpsd(xori,xsft,window,noverlap,nfft,Fs);
%         p(ista)=plot(ax,ft,unwrap(angle(Pxy))/2/pi*Fs./ft,'o-','linewidth', 1,'color',...
%                      color(ista),'markers',1.5);
% %     p(ista)=plot(ax,ft,angle(Pxy)/2/pi*Fs./ft,'o-','linewidth', 1,'color',...
% %                      color(ista),'markers',1.5);
%     plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);    
%     plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
%     plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
% %     ax.XScale = 'log';
%     text(ax,0.05,0.9,'Broadband','unit','normalized','horizontalalignment','left','fontsize',10);
%     text(ax,0.05,0.85,strcat(num2str(sps),{' Hz'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.8,strcat({'nfft: '},num2str(nfft)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.75,strcat({'noverlap: '},num2str(noverlap)),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10);
%     text(ax,0.05,0.7,strcat({'Hann'}),'unit','normalized','horizontalalignment',...
%          'left','fontsize',10); 
%     xlabel(ax,'Frequency (Hz)');
%     ylabel(ax,'Lag (sample)');
% %     lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
%     grid on
%     box on
%     hold(ax, 'off');


% end    

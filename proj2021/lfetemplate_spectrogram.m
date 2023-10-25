%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to examine the spectrogram of the LFE templates, which is essentially
% the time window-averaged PSD. In 'lfetemplate_spectrum.m', we obtained the
% PSD estimate using the 'periodgram' (rather than 'pmtm' or 'pchave'). Since
% the spectrum doesn't tell us the how the power at different frequencies vary
% with the time, especially won't tell us the energy before the dipole, and 
% meanwhile, the bandpass filtering would be numerically unstable when it goes
% to a very narrow passband. So if we really want to know the frequency content
% variation of the other time windows besides the main arrival, we need to obtain
% the spectrogram
%
% -- It seems that spectrogram conveys some information. but not much more than
%   that from the octave-filtering of templates. But surely it is clear from
%   the spectrogram to how low the frequency range that the main arrival can
%   be distinguished from the noise.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/05/05
% Last modified date:   2021/05/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
format short e   % Set the format to 5-digit floating point
clear
clc
% close all

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
    tracebb = stack(:, mid-wlen/2+1: mid+wlen/2)';

    % taper with tukeywin, which is actually a tapered cosine window
    w = tukeywin(size(tracebb,1),0.1);
    for ista = 1:nsta
        tracebb(:,ista) = tracebb(:,ista).*w;        
    end
    
    
    %% obtain the spectrogram
        %%% Recommended this way %%%%%%%%%%%%       
        %%%%%%%%%%%%% the following is to test the function written by FJS
        %%% Well, it seems that 'spectrogram2' basically gives the similar result
        %%% with the same parameter setting.
        seglen = 1024;
        window = hann(seglen);
        nfft = pow2(nextpow2(seglen));
        Fs = sps;
        pover = 0.75;
        noverlap = seglen*pover; 

        octavef = 2.^(3:-1:-5);

        f.fig = figure;
        f.fig.Renderer = 'painters';
        
        widin = 8;  % maximum width allowed is 8.5 inches
        htin = 9;   % maximum height allowed is 11 inches
        set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
        nrow = 3;
        ncol = 1;
        for isub = 1:nrow*ncol
            f.ax(isub) = subplot(nrow,ncol,isub);
            f.ax(isub).Box = 'on';
        end
        color=['r','b','k'];
        
        %%% reposition
        hei = 0.25;
        len = 0.75;
        sep = 0.05;
        for isub = nrow*ncol: -1 :1
            if isub == nrow*ncol
                nhei = 0;
                nsep = 0;
            else
                nhei = nrow*ncol-isub;
                nsep = nrow*ncol-isub;
            end
            set(f.ax(isub), 'position', [ 0.1, 0.08+nhei*hei+nsep*sep, len, hei]);
        end

        for ista = 1: nsta
            [psd, freq, t, psddB]=spectrogram2(tracebb(:,ista),nfft,Fs,seglen,noverlap,'s');
            ax = f.ax(ista);
            hold(ax, 'on');
            
            thr = t+seglen/Fs/2;
%             thr = t;
            surf(ax, thr,freq,10*log10(psd),'EdgeColor','none');
%             imagesc(ax, thr,freq,psddB);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Frequencies bin in log space
%             freqlog = logspace(log10(freq(2)), log10(freq(end)), length(freq) -1);
%             % interpolate the spectrogram
%             psddBi = interp1(freq, psddB, freqlog);
%             imagesc(ax, thr,freqlog,psddBi);
%             ax.YTick = lin2logpos([0.1 1 10], freqlog(1), freqlog(end));
%             ax.YTickLabel = {'0.1'; '1'; '10'};
%             axis xy;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ylim(ax,[min(freq) max(freq)]);
            xlim(ax,[min(thr) max(thr)]);
            ax.YScale = 'log';
            ylabel(ax, 'Frequency (Hz)','fontsize',10);
                    colormap(ax,'jet');
            c=colorbar(ax,'East');
            pos = ax.Position;
            c.Position = [pos(1)+pos(3)+0.02, pos(2), 0.02, pos(4)];
            c.Label.String = 'Spectral density (energy/Hz)';
%             caxis(ax,[-160,0]);
            view(ax,0,90);
            ax.TickDir = 'out';
            text(ax,0.9,0.9,stas(ista,:),'FontSize',12,'unit','normalized');

            for ifreq = 1: length(octavef)

                plot(ax,ax.XLim,[octavef(ifreq) octavef(ifreq)], '--','color',[0.4 0.4 0.4],'linew', 1);
                if octavef(ifreq) < 1
                    fstr = strcat(num2str(round(1/octavef(ifreq))), {' s'});
                else
                    fstr = strcat(num2str(octavef(ifreq)), {' Hz'});
                end
                text(ax,thr(1)+seglen/Fs/4,octavef(ifreq),fstr,'horizontalalignment',...
                    'center','fontsize',8);
                

            end
        
        end
        
        xlabel(f.ax(3),strcat({'Time (s)'}),'fontsize',10);
        titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(Fs),{' Hz, '},{'wlen: '},...
                num2str(seglen),{' nfft: '},num2str(nfft),{', pover: '},num2str(pover),', Hann');
        hdl = title(f.ax(1),titstr,'fontsize',12);
        movev(hdl,0.05);
        
%         %%% save figure
%         print(f.fig,'-dpdf',strcat(temppath,'/BBtempspectrogram_',fam,'_',FLAG,'_',...
%               num2str(sps),'.pdf'));

        
end    
    
    
    
    
    
    
    
    
    
    
    
    
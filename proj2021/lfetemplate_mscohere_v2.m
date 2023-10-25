%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to examine the magnitude-squared coherence, and cross-spectrum of
% the LFE templates. And do some testing to the phase of cross-spectrum
%
% Version 2, difference compared to version 1:
%   the original made templates at different stations are aligned based on 
%   the 4th col of the rotation parameters obtained from 'calc_rots_noresp'
%   which was filtered through 0.5-6.5 hz. So that's why the obtained filtering
%   effect in HF is close to 0 but not exactly 0;
%   Here we use the same template, but as we did in filtering effect correction,
%   we first aligned the templates in LF, then get the phase of the cross-spectrum
%   and see what would happen. (e.g., if there would be 0-slope line through 0,0
%   at LF, meaning the time offset at LF range is constant; and if there would
%   be abrupt boundary at ~1.25 hz; and what would be the HF be like)
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
    
    dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
    flist= ('datalist');
    outf = ('staloc.txt');
    stainfo = getstainfo(dtdir, flist, outf);
    ind=[11,3,4];
    stainfo = stainfo(ind,:);
    
    stalat = str2double(stainfo(:,2));
    stalon = str2double(stainfo(:,3));
    staloc = [stalon stalat zeros(length(stalon),1)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% VERSION 1 for the locations of lfe families, by inversion from hypoinverse
    %%% this is inverted from (0,0) of all fams, same order, location of control points
    %%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
    cont = [
            -123.492667 48.451500 38.1400;
            -123.772167 48.493000 35.5900;
            -123.863167 48.528167 35.2100;
            -123.603333 48.440167 36.7100;
            -123.800167 48.408833 34.5200;
            -123.893333 48.536500 35.0700;
            -123.864500 48.498667 34.8800;
            -123.753333 48.525667 36.2000;
            -123.703667 48.502667 36.4100;
            -123.814333 48.538667 35.7900;
            -123.838500 48.544833 35.6600;
            -123.908000 48.494167 34.5100;       % 006
            -123.879667 48.446167 34.2600;       % 001
            ];
      
    %frequency that has the largest amplitude according to the spectrum, NOT necessarily the best
    %measure, but can be served as a reference
    ftampmaxft = [
                  3.4375e+00   2.4219e+00   4.4531e+00
                  2.1875e+00   2.5000e+00   2.8906e+00
                  2.1875e+00   2.8125e+00   3.2031e+00
                  2.1094e+00   4.0625e+00   4.0625e+00
                  2.1094e+00   2.4219e+00   2.8906e+00
                  2.1094e+00   1.2500e+00   3.2031e+00
                  2.1875e+00   2.8125e+00   3.3594e+00
                  1.9531e+00   1.9531e+00   3.3594e+00
                  3.6719e+00   2.3438e+00   3.9062e+00
                  2.5000e+00   2.7344e+00   3.0078e+01
                  2.1094e+00   2.7344e+00   6.0156e+00
                  2.5000e+00   2.8125e+00   2.1094e+00
                  1.9531e+00   1.5625e+00   1.4844e+00
                  ];
  

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
      
    dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
    flist= ('datalist');
    outf = ('staloc.txt');
    stainfo = getstainfo(dtdir, flist, outf);
    ind=[5,8,6];
    stainfo = stainfo(ind,:);
    
    stalat = str2double(stainfo(:,2));
    stalon = str2double(stainfo(:,3));
    staloc = [stalon stalat zeros(length(stalon),1)];
    
    % location of control fam, this is infered from 'mig_linear_fit_v3' and the detections
    cont = [
            -123.585000 48.436667 36.8800;   % 002   
            -123.549500 48.540833 38.5600;   % 243
            -123.382333 48.574167 40.9800;   % 240
            ];

end
    
% number of used stations
nsta=size(stas,1);         %  number of stations

ALIGN = 'individual';
% ALIGN = 'constrained';
disp(ALIGN)

% PRE = 'lf';
PRE = 'hf';
disp(PRE)


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
    
    
    %% lf filtered and aligned version, similar to the filtering effect part
    switch PRE
        case 'lf'
            wlenseclf = 16;
            wlenlf = wlenseclf*sps;
            lo = 0.5;
            hi = 1.25;
            npo = 2;
            npa = 2;
            cyclskip = 20;
            mshift=14+cyclskip;  % 19/14, samples based on sps in detection, not on 80 sps, but should not affect
            loopoffmax = 4;
            xcmaxAVEnmin = 0.35;  % 0.4/0.35            
            for ista = 1: nsta
                stacktemp(ista,:) = Bandpass(stack(ista,:), sps, lo, hi, npo, npa, 'butter');
            end            
            istart = mid-wlenlf/2;
            iend = istart+wlenlf-1;
            
        case 'hf'
            wlensechf = 4;
            wlenhf = wlensechf*sps;
            lohf = 1.25;
            hihf = 6.5;
            npo = 2;
            npa = 2;
            mshift = 14;    % 29/14, samples based on sps in detection, not on 80 sps, but should not affect
            loopoffmax = 1.5;
            xcmaxAVEnmin = 0.4;
            for ista = 1: nsta
                stacktemp(ista,:) = Bandpass(stack(ista,:), sps, lohf, hihf, npo, npa, 'butter');
            end
            istart = mid-wlenhf/2;
            iend = istart+wlenhf-1;

    end
       
    stackuse = stacktemp(1:3,:);
    stackauto = stackuse.*stackuse;
    lenx = templen-2*mshift;
    stack12x = zeros(lenx, 2*mshift+1);
    stack13x = zeros(lenx, 2*mshift+1);    % stas: 1->PGC, 2->SSIB, 3->SILB
    stack32x = zeros(lenx, 2*mshift+1);
    
    for n=-mshift:mshift
        % PGC corr SSIB, 1+mshift:tracelen-mshift == 1+mshift-n:tracelen-mshift-n == lenx
        stack12x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
            stackuse(2,1+mshift-n:templen-mshift-n);
        % PGC corr SILB
        stack13x(:,n+mshift+1)=stackuse(1,1+mshift:templen-mshift).* ...
            stackuse(3,1+mshift-n:templen-mshift-n);
        % SILB corr SSIB
        stack32x(:,n+mshift+1)=stackuse(3,1+mshift:templen-mshift).* ...
            stackuse(2,1+mshift-n:templen-mshift-n);
    end
    
    sumstack12=zeros(1,2*mshift+1);   % PGC-SSIB
    sumstack13=zeros(1,2*mshift+1);   % PGC-SILB
    sumstack32=zeros(1,2*mshift+1);   % SILB-SSIB
    sumstack1sq=zeros(1,2*mshift+1);  % "sq" are the running cumulative sum, to be used later by differencing the window edpoints
    sumstack2sq=zeros(1,2*mshift+1);
    sumstack3sq=zeros(1,2*mshift+1);
    sumstack3Bsq=zeros(1,2*mshift+1);
    
    sumstack12(1,:) = sum(stack12x(istart-mshift: iend-mshift, :));
    sumstack13(1,:) = sum(stack13x(istart-mshift: iend-mshift, :));
    sumstack32(1,:) = sum(stack32x(istart-mshift: iend-mshift, :));
    sumstack1sq(1,:) = sum(stackauto(1, istart: iend));
    sumstack3Bsq(1,:) = sum(stackauto(3, istart: iend));
    for m = -mshift:mshift
        sumstack2sq(1,m+mshift+1)=sum(stackauto(2, istart-m: iend-m)); %+m??? (yes).
        sumstack3sq(1,m+mshift+1)=sum(stackauto(3, istart-m: iend-m));
    end
    
    %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
    glitches=1.e-7;
    sumstack1sq=max(sumstack1sq,glitches);
    sumstack2sq=max(sumstack2sq,glitches);
    sumstack3sq=max(sumstack3sq,glitches);    % return maximum between A and B
    
    denomstack12n=realsqrt(sumstack1sq.*sumstack2sq);    % Real square root, An error is produced if X is negative
    denomstack13n=realsqrt(sumstack1sq.*sumstack3sq);
    denomstack32n=realsqrt(sumstack3Bsq.*sumstack2sq);
    
    sumstack12n=sumstack12./denomstack12n;   % suffix 'n' means normalized
    sumstack13n=sumstack13./denomstack13n;
    sumstack32n=sumstack32./denomstack32n;
    
    [xcmaxstack12n,imaxstack12]=max(sumstack12n,[],2);   %Integer-offset max cross-correlation
    [xcmaxstack13n,imaxstack13]=max(sumstack13n,[],2);   % along row, max cc val and index in each window
    [xcmaxstack32n,imaxstack32]=max(sumstack32n,[],2);
    
    %Parabolic fit:
    [xmaxstack12n,ymaxstack12n,aSTA12]=parabol(1,mshift,sumstack12n,imaxstack12); %Interpolated max cross-correlation
    [xmaxstack13n,ymaxstack13n,aSTA13]=parabol(1,mshift,sumstack13n,imaxstack13);
    [xmaxstack32n,ymaxstack32n,aSTA32]=parabol(1,mshift,sumstack32n,imaxstack32);
    
    %Center them
    imaxstack12cent=imaxstack12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1: 2*mshift+1, corresponds to -mshift: mshift
    imaxstack13cent=imaxstack13-mshift-1;
    imaxstack32cent=imaxstack32-mshift-1;
    %%% NOTICE: the right order of a closed 3-sta pair is +13, -12, +32, where 13 means 1-->3
    iloopoff=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
    %
    xmaxstack12n=xmaxstack12n-mshift-1;
    xmaxstack13n=xmaxstack13n-mshift-1;
    xmaxstack32n=xmaxstack32n-mshift-1;
    loopoff=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
    loopofffam(ifam) = loopoff;
    xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
    medxcmaxAVEn=median(xcmaxAVEn);
    xmaxstack12ntmp=xmaxstack12n;    % tmp == temporary
    xmaxstack13ntmp=xmaxstack13n;
    xmaxstack32ntmp=xmaxstack32n;
    
    if xcmaxAVEn<xcmaxAVEnmin || abs(loopoff)>loopoffmax || isequal(abs(imaxstack12cent),mshift)...
       || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
        xmaxstack12ntmp=mshift+1; xmaxstack13ntmp=mshift+1; xmaxstack32ntmp=mshift+1; %dummy them, if these criteria are met
        disp('WRONG! The basic criteria is not met');
    else
        xcmaxconprev=-99999.;  %used to be 0; not good with glitches
        %%% NOTE here are using offset imaxSTA12/13/32 without centering, i.e. [1, 2shift+1]
        % Usually, the loop is happened between imaxSTA12n +- floor(loopoffmax+1)
        % width of which is 2*floor(loopoffmax+1)
        %%% floor(2.5)=2; floor(-2.6)=-3
        for iSTA12 = max(1,imaxstack12-floor(loopoffmax+1)): ...
                     min(imaxstack12+floor(loopoffmax+1),2*mshift+1)
            for iSTA13 = max(1,imaxstack13-floor(loopoffmax+1)): ...
                         min(imaxstack13+floor(loopoffmax+1),2*mshift+1)
                ibangon = (mshift+1)-iSTA13+iSTA12;     %%% SEE NOTES #2019/03/17# page 66 to understand better
                %%% i.e., -mshift <= -iSTA13+iSTA12 <= mshift
                if ibangon >= 1 && ibangon <= 2*mshift+1
                    xcmaxcon=sumstack12n(iSTA12)+sumstack13n(iSTA13)+sumstack32n(ibangon);
                    if xcmaxcon > xcmaxconprev
                        xcmaxconprev=xcmaxcon;
                        iSTA12bang=iSTA12;
                        iSTA13bang=iSTA13;
                    end
                end
            end
        end
        %%% will result in the max xcmaxcon and corresponding iSTA12,
        %%% iSTA13, and save them into xcmaxconprev, iSTA12bang and iSTA13bang
        iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
        if abs(iSTA12bang-imaxstack12) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
           abs(iSTA13bang-imaxstack13) <= loopoffmax && ...  % SHOULD not be, i think, could be floor(loopoffmax+1) > loopoffmax
           abs(iSTA32bang-imaxstack32) <= loopoffmax && ...
           sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang) ...
               >= 3*xcmaxAVEnmin   % xcmaxAVEnmin, predetermined
            %%% ALSO, sumsSTA12n(n,iSTA12bang) == sumsSTA12nn(iSTA12bang)
            
            xmaxstack12ntmp=iSTA12bang-(mshift+1); %without interpolation this is just centering.
            xmaxstack13ntmp=iSTA13bang-(mshift+1);
            xmaxstack32ntmp=iSTA32bang-(mshift+1);
            
            %%% xcmaxAVEnbang is added by Chao, to distinguish from
            %%% xcmaxAVEn, because it is max average CC coef
            xcmaxAVEnbang=(sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang))/3;
        end
    end
    
    if xmaxstack13ntmp-xmaxstack12ntmp+xmaxstack32ntmp ~=0
        disp('WRONG! Loopoff is not enclosed');
    else
%         disp('Loopoff in lf is 0');
    end
    
    
    %define the new desired length of templates, above length is for obtaining the alignment only
    wlensec = 16;
    wlen = wlensec*sps;   
    istart = mid-wlen/2;
    iend = istart+wlen-1;

    % align the broadband stacked records based on alignment in LF
    switch ALIGN
        case 'individual'
            x(:, 1) = stack(1, istart: iend);
            x(:, 2) = stack(2, istart-imaxstack12cent: iend-imaxstack12cent);
            x(:, 3) = stack(3, istart-imaxstack13cent: iend-imaxstack13cent);
        case 'constrained'
            x(:, 1) = stack(1, istart: iend);
            x(:, 2) = stack(2, istart-xmaxstack12ntmp: iend-xmaxstack12ntmp);
            x(:, 3) = stack(3, istart-xmaxstack13ntmp: iend-xmaxstack13ntmp);
    end
    
    
    
    %% phase of cross spectrum in normalized phase lag
%     nfft = pow2(nextpow2(size(x,1))-1);
    nfft = min(2048, pow2(nextpow2(size(x,1))-1));
    window = hann(nfft);
    noverlap = nfft*2/4;
    Fs = sps;
    % filter param
    lolf = 0.5;
    lohf = 1.25;
    hihf = 6.5;
    
    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 7.5;  % maximum width allowed is 8.5 inches
    htin = 9;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 3;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        f.ax(isub).FontSize = 9;
        f.ax(isub).Box = 'on';
        grid(f.ax(isub), 'on');
        f.ax(isub).GridLineStyle = '--';
    end
        
    % reposition
    set(f.ax(3), 'position', [ 0.1, 0.1, 0.85, 0.25]);
    set(f.ax(2), 'position', [ 0.1, 0.4, 0.85, 0.25]);
    set(f.ax(1), 'position', [ 0.1, 0.7, 0.85, 0.25]);

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
        end
%         ylim(ax,[-2.5,0.5]);
%         for ista = 1: nsta
% %             if ista == 2
%                 p(ista)=plot(ax,ft,unwrap(angle(Pxy(:,ista)))/pi,'o-','linewidth', 1,'color',...
%                              color(ista),'markers',1.5);
% %             else
% %                 p(ista)=plot(ax,ft,angle(Pxy(:,ista))/pi,'o-','linewidth', 1,'color',...
% %                              color(ista),'markers',1.5); 
% %             end
%         end        
        plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        xlabel(ax,'Frequency (Hz)', 'fontsize',11);
        ylabel(ax,'Normalized Phase lag (rad)','fontsize',11);
        lgd = legend(ax,p,{'TWKB-LZB','TWKB-MGCB','LZB-MGCB'},'location','southwest','fontsize',8);
        grid(ax, 'on');
        ax.Box = 'on';
        hold(ax, 'off');
        
    end
    titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},num2str(wlensec),...
                    {' s, '},{'nfft: '},num2str(nfft),{', noverlap: '},num2str(noverlap),...
                    ', Hann');
    hdl = title(f.ax(1),titstr,'fontsize',12);
    movev(hdl,0.05);

%     %%% save figure
%     print(f.fig,'-dpdf',strcat(temppath,'/BBtempCPSD_',fam,'_',FLAG,'_',PRE,'align_',...
%           num2str(sps),'_',num2str(wlensec),'s.pdf'));
    close(f.fig)

    %% have a comparison with the prediction from a freq-independent attenuation
    % lines of codes are from 'testattenuationmodel.m'
    [starel(:,1),starel(:,2)] = absloc2relaloc(stalon,stalat,cont(ifam,1),cont(ifam,2));
    starel(:,3) = ones(size(starel,1),1)*cont(ifam,3);
    Dist = sqrt(sum(starel.^2, 2));  % distance from source to any single station

    %let's say the velocity model gives us the vel. at 1 Hz
    c0 = 3.84;
    f0 = 1;
    om0 = 2*pi*f0;
    Fs = sps;
    Qs = 20; 
    
    %if fix a frequency, and analyse the tshift/nshift, or normalized phase as a function of 
    %frequency, meanning that you align the records at the reference frequency at two stations, and
    %then analyze the time offset at other frequency wrt the aligned frequency 
    falign = [(1.25+6.5)/2 2.5];
%     falign = 2.5;
    freq = 0.01: 0.02: 8;

    f.fig = figure;
    f.fig.Renderer = 'painters';
    widin = 7.5;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);
    nrow = 3;
    ncol = 1;
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        f.ax(isub).FontSize = 9;
        f.ax(isub).Box = 'on';
        grid(f.ax(isub), 'on');
        f.ax(isub).GridLineStyle = '--';
    end
        
    % reposition
    set(f.ax(3), 'position', [ 0.1, 0.1, 0.85, 0.25]);
    set(f.ax(2), 'position', [ 0.1, 0.4, 0.85, 0.25]);
    set(f.ax(1), 'position', [ 0.1, 0.7, 0.85, 0.25]);
    
    color = ['r';'b';'k'];
%     color = ['k';'k';'k'];
    stastr = {'TWKB-LZB';'TWKB-MGCB';'LZB-MGCB'};
    lines = {'--'; '-.'};
    linew = [1.5 1];

    for ista = 1: nsta
        ax = f.ax(ista);
        hold(ax, 'on');
        xlim(ax,[0, 8]);
        ylim(ax,[-0.6,0.6]);
        yticks(ax, -0.6: 0.2: 0.6);
        plot(ax,[lolf lolf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        plot(ax,[lohf lohf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        plot(ax,[hihf hihf],ax.YLim, '--','color',[0.6 0.6 0.6],'linew', 1.5);
        p(1)=plot(ax,ft,angle(Pxy(:,ista))/pi,'o-','linewidth', 1,'color',...
            color(ista),'markers',2);
        text(ax,0.2,0.9,stastr{ista},'fontsize',11,'unit','normalized',...
            'Horizontalalignment','left','EdgeColor','k','Margin',2);
        for i = 1: length(falign)
            omalign = 2*pi*falign(i);
            omlf = 2*pi*freq;
            calign = c0*(1+log(omalign/om0)/(pi*Qs));
            dtfreq_12 = (Dist(2) - Dist(1))*log(omalign./omlf)./(pi*Qs*calign);
            dtfreq_13 = (Dist(3) - Dist(1))*log(omalign./omlf)./(pi*Qs*calign);
            dtfreq_23 = (Dist(3) - Dist(2))*log(omalign./omlf)./(pi*Qs*calign);

            phalagfreq_12 = dtfreq_12 .* omlf /pi;
            phalagfreq_13 = dtfreq_13 .* omlf /pi;
            phalagfreq_23 = dtfreq_23 .* omlf /pi;

            if ista == 1
                p(i+1) = plot(ax,freq,phalagfreq_12,lines{i},'linewidth',linew(i),'color',...
                                color(ista));
            elseif ista == 2
                p(i+1) = plot(ax,freq,phalagfreq_13,lines{i},'linewidth',linew(i),'color',...
                                color(ista));
            elseif ista == 3
                p(i+1) = plot(ax,freq,phalagfreq_23,lines{i},'linewidth',linew(i),'color',...
                                color(ista));
            end

        end
        lgd = legend(ax,p,{'Obs. by aligning in HF passband',['Pred. by aligning at ',...
            num2str(falign(1)),' Hz'], ['Pred. by aligning at ',...
            num2str(falign(2)),' Hz']},'fontsize',8);
        
        if ifam == 1
            if ista == 1
                lgd.Location = 'south';
            elseif ista == 2
                lgd.Location = 'southwest';
            else
                lgd.Location = 'south';
            end
        elseif ifam == 6
            lgd.Location = 'southwest';

        end
                        
        % make background transparent
        set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

        ylabel(ax,'Normalized Phase lag (rad)','fontsize', 11);
    end
    xlabel(f.ax(3),'Frequency (Hz)','fontsize', 11);
    
    text(f.ax(1),0.5,0.9,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
    text(f.ax(1),0.5,0.75,strcat({'Q_s: '},num2str(Qs)),'fontsize',11,'unit','normalized',...
         'Horizontalalignment','left');

%     titstr = strcat(FLAG,{', '},fam,{', Broadband, '},num2str(sps),{' Hz, '},num2str(wlensec),...
%                     {' s, '},{'nfft: '},num2str(nfft),{', noverlap: '},num2str(noverlap),...
%                     {', Hann, '},{'Qs: '},num2str(Qs));
%     hdl = title(f.ax(1),titstr,'fontsize',12);
%     movev(hdl,0.05);
    
    
    %%% save figure
    print(f.fig,'-dpdf',strcat(temppath,'/BBtempCPSD_vs_atten_',fam,'_',FLAG,'_',PRE,'align_',...
          num2str(sps),'_',num2str(wlensec),'s.pdf'));
%     close(f.fig)

end

  











  

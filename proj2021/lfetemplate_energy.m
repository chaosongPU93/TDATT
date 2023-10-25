%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to examine the energy of the main arrival, energy ratio between the
% main arival and surrounding windows, and hilbert envelope of the LFE templates
% The purpose behind is to find some thresholds for the isolated LF arrivals.
% 'isolate' can means the energy ratio between the strongest arrival
% window and its prior window and posterior window.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/03/25
% Last modified date:   2021/03/25
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
sps = 40;
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
    
    
    %% plot the templates
    %%% plot the 1-step template
    %%% figure 3
    wlensec = 16;
    wlen = wlensec*sps;
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    for ista = 1: nsta
        plot(ax,(1:wlen)/sps,stack(ista, mid-wlen/2+1: mid+wlen/2), 'linewidth', 1,'color',...
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
    xlim(ax,[0 wlensec]);
%     ylim(ax,[0 4]);
    ax.Box='on';
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    
    x = stack(:, mid-wlen/2+1: mid+wlen/2)';
    
    %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
    [envup, envlo] = envelope(x);
    
    medenvup = median(envup,2);
    medenvlo = median(envlo,2);
    
    plot(ax,(1:wlen)/sps,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
    plot(ax,(1:wlen)/sps,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
    hold(ax, 'off');
    
            
    %% lf filtered and aligned version, similar to the filtering effect part
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
    
    istart = mid-wlenlf/2;
    iend = istart+wlenlf-1;
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
    iloopofflf=imaxstack13cent-imaxstack12cent+imaxstack32cent; %How well does the integer loop close?
    %
    xmaxstack12n=xmaxstack12n-mshift-1;
    xmaxstack13n=xmaxstack13n-mshift-1;
    xmaxstack32n=xmaxstack32n-mshift-1;
    loopoff=xmaxstack13n-xmaxstack12n+xmaxstack32n; %How well does the interpolated loop close?
    loopofflf(ifam) = loopoff;
    xcmaxAVEn=(xcmaxstack12n+xcmaxstack13n+xcmaxstack32n)/3;
    medxcmaxAVEn=median(xcmaxAVEn);
    xmaxstack12ntmplf=xmaxstack12n;    % tmp == temporary
    xmaxstack13ntmplf=xmaxstack13n;
    xmaxstack32ntmplf=xmaxstack32n;
    
    if xcmaxAVEn<xcmaxAVEnmin || abs(loopoff)>loopoffmax || isequal(abs(imaxstack12cent),mshift)...
       || isequal(abs(imaxstack13cent),mshift) || isequal(abs(imaxstack32cent),mshift)
        xmaxstack12ntmplf=mshift+1; xmaxstack13ntmplf=mshift+1; xmaxstack32ntmplf=mshift+1; %dummy them, if these criteria are met
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
            
            xmaxstack12ntmplf=iSTA12bang-(mshift+1); %without interpolation this is just centering.
            xmaxstack13ntmplf=iSTA13bang-(mshift+1);
            xmaxstack32ntmplf=iSTA32bang-(mshift+1);
            
            %%% xcmaxAVEnbang is added by Chao, to distinguish from
            %%% xcmaxAVEn, because it is max average CC coef
            xcmaxAVEnbang=(sumstack12n(iSTA12bang)+sumstack13n(iSTA13bang)+sumstack32n(iSTA32bang))/3;
        end
    end
    
    if xmaxstack13ntmplf-xmaxstack12ntmplf+xmaxstack32ntmplf ~=0
        disp('WRONG! Loopoff is not enclosed');
    else
%         disp('Loopoff in lf is 0');
    end
    
    % align the records
    xlf(:, 1) = stackuse(1, istart: iend);
    xlf(:, 2) = stackuse(2, istart-xmaxstack12ntmplf: iend-xmaxstack12ntmplf);
    xlf(:, 3) = stackuse(3, istart-xmaxstack13ntmplf: iend-xmaxstack13ntmplf);
    
    % find the arrival time of the strongest concentration window
    STA12tr=xlf(:, 1).* xlf(:, 2); 
    STA13tr=xlf(:, 1).* xlf(:, 3);
    STA32tr=xlf(:, 3).* xlf(:, 2);
    cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr);    % sum of the cumsum of all traces
    %%% first get the squared sum of each 0.5s window, then get the maximum and start indice
    %%% this gives the proxy of Energy of 0.5s window, i.e. the integral of \dot(E)
    %%% of eqn 1 in Rubin et al. 2013
    %%% idiff<==>win [idiff+1, idiff+20] in regional ind<==>win [istart+idiff, istart+idiff+19] in global index
    concentration=1.1; %in seconds; how concentrated is the coherent energy within the window?
    cncntr=concentration*sps;   % in samples, 20
    [cumsumtrdiff, idiff]=max(cumsumtr(cncntr+1:wlenlf)-cumsumtr(1:wlenlf-cncntr));

    
    %%
    f.fig = figure;
    f.fig.Renderer = 'painters';
    color=['r','b','k'];
    ax = gca;
    hold(ax, 'on');
    
    maxamppos = max(xlf, [],1);
    minampneg = min(xlf, [],1);
    yma=max(maxamppos);
    ymi=min(minampneg);
    yma=2*max(yma,-ymi);

    %%%%%%%% normalize each trace by its own maximum
    maxamp = max([maxamppos; -minampneg], [], 1);
    for ii = 1: 3
        xlf(:, ii) = xlf(:, ii)./maxamp(ii);
    end
    yma=2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    xarea=[idiff/sps (idiff+cncntr)/sps (idiff+cncntr)/sps idiff/sps idiff/sps];
    yarea=[-yma -yma yma yma -yma];
    patch(ax,xarea,yarea,'k','Facealpha',0.1,'edgecolor','none');
  
    for ista = 1: nsta
        plot(ax,(1:wlenlf)/sps, xlf(:, ista), 'linewidth', 1,'color',color(ista));
%         text(ax,1,0,stas(ista,:),'fontsize',9);
    end
    text(ax,0.95,0.9,'Aligned in 0.5-1.25 Hz','fontsize',10,'unit','normalized','Horizontalalignment','right');
    text(ax,0.03,0.92,'a','FontSize',10,'unit','normalized','EdgeColor','k','Margin',2);
    text(ax,0.03,0.1,fam,'fontsize',12,'unit','normalized','Horizontalalignment','left');
    text(ax,0.95,0.1,FLAG,'fontsize',12,'unit','normalized','Horizontalalignment','right');
    
    xlabel(ax,'Time (s)','fontsize',10);
%     xlabel(ax,strcat('Samples (at',{' '},num2str(sps),{' '}, 'Hz)'),'fontsize',10);
    ylabel(ax,'Amplitude','fontsize',10);
    xlim(ax,[0 wlenseclf]);
    ylim(ax,[-yma yma]);
    ax.Box='on';
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    
    %%%%%%% Below: Find the max of the hilbert transform within the arrival; then walk out from there +/-
    %     wavehil = abs(hilbert(tracefile(winlen*(ndet-1)+1:winlen*ndet,2:4)));
    [envuplf, envlolf] = envelope(xlf);
    
    medenvup = median(envuplf,2);
    medenvlo = median(envlolf,2);
    
%     plot(ax,time(winlen*(ndet-1)+1:winlen*ndet),medhil,'color',[0.3 0.3 0.3],'linew',0.6); %just to test Hilbert transform
    plot(ax,(1:wlenlf)/sps,medenvup,'color',[0.7 0.7 0.7],'linew',1.5);
    plot(ax,(1:wlenlf)/sps,medenvlo,'color',[0.7 0.7 0.7],'linew',1.5);
    
    [maxhil, loc]=max(medenvup(idiff: idiff+cncntr-1));
    maxloc=idiff+loc-1;
    icheck=maxloc;
    while medenvup(icheck)>0.5*maxhil && icheck<wlenlf-1
        icheck=icheck+1;
    end
    hilwid(ifam,1)=icheck-maxloc;
    while medenvup(icheck)>0.25*maxhil && icheck<wlenlf-1
        icheck=icheck+1;
    end
    hilwid(ifam,2)=icheck-maxloc;
    while medenvup(icheck)>0.125*maxhil && icheck<wlenlf-1
        icheck=icheck+1;
    end
    hilwid(ifam,3)=icheck-maxloc;

    icheck=maxloc;
    while medenvup(icheck)>0.5*maxhil && icheck>1
        icheck=icheck-1;
    end
    hilwid(ifam,4)=maxloc-icheck;
    while medenvup(icheck)>0.25*maxhil && icheck>1
        icheck=icheck-1;
    end
    hilwid(ifam,5)=maxloc-icheck;
    while medenvup(icheck)>0.125*maxhil && icheck>1
        icheck=icheck-1;
    end
    hilwid(ifam,6)=maxloc-icheck;
    
    % This for longer-term averages:
    maxhil=median(medenvup(idiff: idiff+cncntr-1)); %median of arrival window
    maxloc=idiff; %"maxloc" is really start of arrival window
    icheck=maxloc;
    while median(medenvup(icheck:icheck+cncntr-1))>0.666*maxhil && icheck+cncntr<wlenlf 
        icheck=icheck+1;
    end
    hilwid(ifam,7)=icheck-maxloc;
    while median(medenvup(icheck:icheck+cncntr-1))>0.5*maxhil && icheck+cncntr<wlenlf 
        icheck=icheck+1;
    end
    hilwid(ifam,8)=icheck-maxloc;
    while median(medenvup(icheck:icheck+cncntr-1))>0.25*maxhil && icheck+cncntr<wlenlf 
        icheck=icheck+1;
    end
    hilwid(ifam,9)=icheck-maxloc;
    
    maxloc=idiff+cncntr-1; %"maxloc" is really end of arrival window
    icheck=maxloc;
    while median(medenvup(icheck-cncntr+1:icheck))>0.666*maxhil && icheck-cncntr>0 
        icheck=icheck-1;
    end
    hilwid(ifam,10)=maxloc-icheck;
    while median(medenvup(icheck-cncntr+1:icheck))>0.5*maxhil && icheck-cncntr>0 
        icheck=icheck-1;
    end
    hilwid(ifam,11)=maxloc-icheck;
    while median(medenvup(icheck-cncntr+1:icheck))>0.25*maxhil && icheck-cncntr>0
        icheck=icheck-1;
    end
    hilwid(ifam,12)=maxloc-icheck;
    
    plot(ax.XLim, [0.5*maxhil 0.5*maxhil], '--');
    plot(ax.XLim, [0.25*maxhil 0.25*maxhil], '--');
    plot(ax.XLim, [0.125*maxhil 0.125*maxhil], '--');
    
    hold(ax, 'off');

    %%%%%%% Above: Find the max of the hilbert transform within the arrival; then walk out from there +/-

    %%        
    %%%%%%% Below: Find the global and nearest local maximum and minimum, compare the amplitude
    maxamppos = max(xlf, [],1);
    minampneg = min(xlf, [],1);
    for ii = 1: 3
        gmaxloc = find(xlf(:, ii) == maxamppos(ii));
        [lmaxamp(ii,4:5),loc] = findpeaks(xlf(gmaxloc+1: end, ii),'NPeaks',2,...
            'MinPeakDistance', 2);
        lmaxloc(ii,4:5) = loc+gmaxloc;
        [lmaxamp(ii,2:-1:1),loc] = findpeaks(flipud(xlf(1: gmaxloc-1, ii)),'NPeaks',2,...
            'MinPeakDistance', 2);
        lmaxloc(ii,2:-1:1) = gmaxloc-loc;
        lmaxamp(ii,3) = maxamppos(ii);
        lmaxloc(ii,3) = gmaxloc;
        
        gminloc = find(xlf(:, ii) == minampneg(ii));
        [pk,loc] = findpeaks(-xlf(gminloc+1: end, ii),'NPeaks',2,...
            'MinPeakDistance', 2);
        lminamp(ii,4:5) = -pk;
        lminloc(ii,4:5) = loc+gminloc;
        [pk,loc] = findpeaks(flipud(-xlf(1: gminloc-1, ii)),'NPeaks',2,...
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
    medlmaxarat(ifam, :) = median(lmaxarat, 1);
    medlminarat(ifam, :) = median(lminarat, 1);
    lmaxwid(:, 1:2) = lmaxloc(:,3)-lmaxloc(:, 1:2);
    lmaxwid(:, 3:4) = lmaxloc(:,4:5)-lmaxloc(:, 3);
    lminwid(:, 1:2) = lminloc(:,3)-lminloc(:, 1:2);
    lminwid(:, 3:4) = lminloc(:,4:5)-lminloc(:, 3);
    medlmaxwid(ifam, :) = median(lmaxwid, 1);
    medlminwid(ifam, :) = median(lminwid, 1);
    
    %%%%%%% Above: Find the global and nearest local maximum and minimum, compare the amplitude
        
    
    %%
    %%%%%%%%%%% below is to find the energy of different windows	%%%%%%%%%%%%%%%%%%%%%%
    %What is amp squared in strongest coherent 1/2 sec?
    %%% amp squared is a self to self operation
    %%% NOTICE! HERE, istart=igstart+(n-1)*winoff >= igstart,
    %%% So, isdiff >= istart >= igstart
    offset=round(0.5*cncntr);   % +- 1/2*cncntr, 10 samples, 0.25s
    isdiff=idiff; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
    iediff=idiff-1+cncntr;  % cncntr is 20sps == 0.5s
        
    dummy=xlf(isdiff: iediff, 1).^2+ ...   % point square
          xlf(isdiff: iediff, 2).^2+ ...
          xlf(isdiff: iediff, 3).^2;
    dum2=cumsum(dummy);
    %                     Ampsq(nin+1)=dum2(end);   % Ampsq(1) == amplitude square = cumsum of all dummy = squared sum
    Ampsq=dum2(end) / length(dummy);
    
    %%Energy in prior 2.5*cncntr seconds, with offset (assuming 0.5cncntr)
    %%% energy is a self to self operation, and is proportional to amp squared
    %%% the offset could be regarded as a buffer zone
    if isdiff > round(2.5*cncntr)+(mshift-cyclskip)+offset
        dummy=xlf(isdiff-round(2.5*cncntr)-offset: isdiff-offset-1, 1).^2+ ...
              xlf(isdiff-round(2.5*cncntr)-offset: isdiff-offset-1, 2).^2+ ...
              xlf(isdiff-round(2.5*cncntr)-offset: isdiff-offset-1, 3).^2;
    else
        dummy=xlf((mshift-cyclskip): isdiff-1-offset, 1).^2+ ...
              xlf((mshift-cyclskip): isdiff-1-offset, 2).^2+ ...
              xlf((mshift-cyclskip): isdiff-1-offset, 3).^2;
    end
    dum2=cumsum(dummy);
    %                     Prev(nin+1)=dum2(end);    % Prev(1) == previous 1.25s window before the strongest window, length 2.5*cncntr
    Prior=dum2(end) / length(dummy);
    clear dummy
       
    %Energy in prior 4 s
    if isdiff > 4*sps +(mshift-cyclskip)+offset
        dummy=xlf(isdiff-4*sps-offset: isdiff-offset-1, 1).^2+ ...
              xlf(isdiff-4*sps-offset: isdiff-offset-1, 2).^2+ ...
              xlf(isdiff-4*sps-offset: isdiff-offset-1, 3).^2;
    else
        dummy=xlf((mshift-cyclskip): isdiff-1-offset, 1).^2+ ...
              xlf((mshift-cyclskip): isdiff-1-offset, 2).^2+ ...
              xlf((mshift-cyclskip): isdiff-1-offset, 3).^2;
    end
    dum2=cumsum(dummy);
    %                     Prior4(nin+1)=dum2(end);
    Prior4=dum2(end) / length(dummy);
    clear dummy
        
    % energy (amp squared sum) in post 2.5*cncntr second window, with offset
    if iediff+round(2.5*cncntr)+(mshift-cyclskip)+offset <= size(xlf,1)
        dummy=xlf(iediff+1+offset: iediff+round(2.5*cncntr)+offset, 1).^2+ ...
              xlf(iediff+1+offset: iediff+round(2.5*cncntr)+offset, 2).^2+ ...
              xlf(iediff+1+offset: iediff+round(2.5*cncntr)+offset, 3).^2;
    else
        dummy=xlf(iediff+1+offset: size(xlf,1)-(mshift-cyclskip), 1).^2+ ...
              xlf(iediff+1+offset: size(xlf,1)-(mshift-cyclskip), 2).^2+ ...
              xlf(iediff+1+offset: size(xlf,1)-(mshift-cyclskip), 3).^2;
    end
    dum2=cumsum(dummy);
    %                     Post(nin+1)=dum2(end);
    Post=dum2(end) / length(dummy);
    clear dummy
    
    % energy (amp squared sum) in post 4 s window, with offset
    if iediff+4*sps+(mshift-cyclskip)+offset <= size(xlf,1)
        dummy=xlf(iediff+1+offset: iediff+4*sps+offset, 1).^2+ ...
              xlf(iediff+1+offset: iediff+4*sps+offset, 2).^2+ ...
              xlf(iediff+1+offset: iediff+4*sps+offset, 3).^2;
    else
        dummy=xlf(iediff+1+offset: size(xlf,1)-(mshift-cyclskip), 1).^2+ ...
              xlf(iediff+1+offset: size(xlf,1)-(mshift-cyclskip), 2).^2+ ...
              xlf(iediff+1+offset: size(xlf,1)-(mshift-cyclskip), 3).^2;
    end
    dum2=cumsum(dummy);
    %                     Post4(nin+1)=dum2(end);
    Post4=dum2(end) / length(dummy);
    clear dummy
    
    
    %Energy in prior half of winlen, i.e. 8 s in current setting
    istart = mid-wlenlf/2;
    iend = istart+wlenlf-1;
    tmplf(:, 1) = stackuse(1, istart: iend);
    tmplf(:, 2) = stackuse(2, istart-xmaxstack12ntmplf: iend-xmaxstack12ntmplf);
    tmplf(:, 3) = stackuse(3, istart-xmaxstack13ntmplf: iend-xmaxstack13ntmplf);
    STA12tr=tmplf(:, 1).* tmplf(:, 2); 
    STA13tr=tmplf(:, 1).* tmplf(:, 3);
    STA32tr=tmplf(:, 3).* tmplf(:, 2);
    cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr);    % sum of the cumsum of all traces
    %%% first get the squared sum of each 0.5s window, then get the maximum and start indice
    %%% this gives the proxy of Energy of 0.5s window, i.e. the integral of \dot(E)
    %%% of eqn 1 in Rubin et al. 2013
    %%% idiff<==>win [idiff+1, idiff+20] in regional ind<==>win [istart+idiff, istart+idiff+19] in global index
    concentration=1.1; %in seconds; how concentrated is the coherent energy within the window?
    cncntr=concentration*sps;   % in samples, 20
    [~, tmpidiff]=max(cumsumtr(cncntr+1:wlenlf)-cumsumtr(1:wlenlf-cncntr));

    tmpisdiff=tmpidiff; %Start of strongest 0.5s, DELETE -1 by Chao, 2019/02/17, see NOTES
    tmpiediff=tmpidiff-1+cncntr;  % cncntr is 20sps == 0.5s
    
    if tmpisdiff > wlenlf/2 +(mshift-cyclskip)+offset
        dummy=tmplf(tmpisdiff-wlenlf/2-offset: tmpisdiff-offset-1, 1).^2+ ...
              tmplf(tmpisdiff-wlenlf/2-offset: tmpisdiff-offset-1, 2).^2+ ...
              tmplf(tmpisdiff-wlenlf/2-offset: tmpisdiff-offset-1, 3).^2;
    else
        dummy=tmplf((mshift-cyclskip): tmpisdiff-1-offset, 1).^2+ ...
              tmplf((mshift-cyclskip): tmpisdiff-1-offset, 2).^2+ ...
              tmplf((mshift-cyclskip): tmpisdiff-1-offset, 3).^2;
    end
    dum2=cumsum(dummy);
    %                     Prior8(nin+1)=dum2(end);
    Prior8=dum2(end) / length(dummy);
    clear dummy

    % energy (amp squared sum) in post winlen/2 window, i.e. 8 s, with offset
    if tmpiediff+wlenlf/2+(mshift-cyclskip)+offset <= size(tmplf,1)
        dummy=tmplf(tmpiediff+1+offset: tmpiediff+wlenlf/2+offset, 1).^2+ ...
              tmplf(tmpiediff+1+offset: tmpiediff+wlenlf/2+offset, 2).^2+ ...
              tmplf(tmpiediff+1+offset: tmpiediff+wlenlf/2+offset, 3).^2;
    else
        dummy=tmplf(tmpiediff+1+offset: size(tmplf,1)-(mshift-cyclskip), 1).^2+ ...
              tmplf(tmpiediff+1+offset: size(tmplf,1)-(mshift-cyclskip), 2).^2+ ...
              tmplf(tmpiediff+1+offset: size(tmplf,1)-(mshift-cyclskip), 3).^2;
    end
    dum2=cumsum(dummy);
    %                     Post8(nin+1)=dum2(end);
    Post8=dum2(end) / length(dummy);
    clear dummy
    
    energy(ifam,:) = [Ampsq Prior Post Prior4 Post4 Prior8 Post8];
    
    eratprior = Ampsq./Prior;
    eratprior4 = Ampsq./Prior4;
    eratprior8 = Ampsq./Prior8;
    eratpost = Ampsq./Post;
    eratpost4 = Ampsq./Post4;
    eratpost8 = Ampsq./Post8;
    erat(ifam, :) = [eratprior eratpost eratprior4 eratpost4 eratprior8 eratpost8];

    %%%%%%%%%%% above is to find the energy of different windows	%%%%%%%%%%%%%%%%%%%%%%
                    
                    
end    


%%
if strcmp(FLAG, 'TWKB') 
    % families in the NW region
    famintind = [2,3,6,7,8,9,10,11,12];
    % disp(nfampool(famintind, :));
    medhilwid = median(hilwid(famintind, :), 1);
    medenergy = median(energy(famintind, :), 1);
    mederat = median(erat(famintind, :), 1);
    medarat = median([medlmaxarat(famintind, :) medlminarat(famintind, :)], 1);
elseif strcmp(FLAG, 'PGC')    
    medhilwid = median(hilwid, 1);
    medenergy = median(energy, 1);
    mederat = median(erat, 1);
    medarat = median([medlmaxarat medlminarat], 1);
end    
hilwidsave = [hilwid; medhilwid];
energysave = [energy; medenergy];
eratsave = [erat; mederat];
aratsave = [medlmaxarat medlminarat; medarat];
    
% write into file    
fid = fopen(strcat(temppath, '/temp_hilbwidth.',FLAG,'.sps',num2str(sps),'.wlen',num2str(wlenseclf), ...
            '.f',num2str(lo),'-',num2str(hi)),'w+');
fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d \n', hilwidsave');
fclose(fid);
    
fid = fopen(strcat(temppath, '/temp_energy.',FLAG,'.sps',num2str(sps),'.wlen',num2str(wlenseclf), ...
            '.f',num2str(lo),'-',num2str(hi)),'w+');
fprintf(fid,'%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n', energysave');
fclose(fid);    
    
fid = fopen(strcat(temppath, '/temp_engratio.',FLAG,'.sps',num2str(sps),'.wlen',num2str(wlenseclf), ...
            '.f',num2str(lo),'-',num2str(hi)),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f \n', eratsave');
fclose(fid);

fid = fopen(strcat(temppath, '/temp_ampratio.',FLAG,'.sps',num2str(sps),'.wlen',num2str(wlenseclf), ...
            '.f',num2str(lo),'-',num2str(hi)),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n', aratsave');
fclose(fid);
    
    
    
    
    
    
    
    
    
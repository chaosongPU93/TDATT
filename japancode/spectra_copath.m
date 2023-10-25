function [pmeq1EssA,pmeq1NssA,pmeq1UssA,pmeq1EpsA,pmeq1NpsA,pmeq1UpsA,pmft,f1,f2,f3]= ...
    spectra_copath(eq1EpsA,eq1NpsA,eq1UpsA,eq1EssA,eq1NssA,eq1UssA,eq1EnoiA,eq1NnoiA,eq1UnoiA,...
    sps,tmpevt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is to read in the P and S signal and noise window data of E
% N and U component data of hinet for all co-path events tied to each target
% single target co-path-time event
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/26
% Last modified date:   2020/07/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);


%% obtain the spectra of different data
% For PSD, use pmtm or periodgram (no averaging) for regular earthquakes; use pchave for tremor and
% noise;
% For amplitude spectrum, use fft, see https://www.mathworks.com/help/matlab/ref/fft.html

%%%% PSD estimate with periodgram
%%% EQ 1
% periodogram, P and S signal

%for i = size(eq1Ess,2)

% nfft = pow2(nextpow2(size(eq1EssA,1))-1);
nfft = pow2(nextpow2(size(eq1EssA,1)));
window = hann(size(eq1EssA,1));
Fs = sps;
[pdeq1EssA,pdft] = periodogram(eq1EssA,window,nfft,Fs);  % psd is done to each col, i.e., each station
[pdeq1NssA,~] = periodogram(eq1NssA,window,nfft,Fs);
[pdeq1UssA,~] = periodogram(eq1UssA,window,nfft,Fs);

[pdeq1EpsA,~] = periodogram(eq1EpsA,window,nfft,Fs);  % psd is done to each col, i.e., each station
[pdeq1NpsA,~] = periodogram(eq1NpsA,window,nfft,Fs);
[pdeq1UpsA,~] = periodogram(eq1UpsA,window,nfft,Fs);

% pmtm, P and S signal
nw = 4;
[pmeq1EssA,pmft] = pmtm(eq1EssA,nw,nfft,Fs);
[pmeq1NssA,~] = pmtm(eq1NssA,nw,nfft,Fs);
[pmeq1UssA,~] = pmtm(eq1UssA,nw,nfft,Fs);

[pmeq1EpsA,~] = pmtm(eq1EpsA,nw,nfft,Fs);
[pmeq1NpsA,~] = pmtm(eq1NpsA,nw,nfft,Fs);
[pmeq1UpsA,~] = pmtm(eq1UpsA,nw,nfft,Fs);


%%% EQ 1, noise around origin
% periodogram
% nfft = pow2(nextpow2(size(eq1EnoiA,1))-1);
nfft = pow2(nextpow2(size(eq1EnoiA,1)));
window = hann(size(eq1EnoiA,1));
Fs = sps;
[pdeq1EnoiA,~] = periodogram(eq1EnoiA,window,nfft,Fs);
[pdeq1NnoiA,~] = periodogram(eq1NnoiA,window,nfft,Fs);
[pdeq1UnoiA,~] = periodogram(eq1UnoiA,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq1EnoiA,~] = pmtm(eq1EnoiA,nw,nfft,Fs);
[pmeq1NnoiA,~] = pmtm(eq1NnoiA,nw,nfft,Fs);
[pmeq1UnoiA,~] = pmtm(eq1UnoiA,nw,nfft,Fs);


%% distinguish those events based on their magnitudes
% cause larger mag events likely have a higher absolute
% PSD
mag = tmpevt(:,11);
ind1 = find(mag<-1);
ind10 = find(mag>=-1 & mag<0);
ind01 = find(mag>=0 & mag<1);
ind12 = find(mag>=1 & mag<2);
ind23 = find(mag>=2 & mag<3);
ind34 = find(mag>=3 & mag<4);
ind45 = find(mag>=4 & mag<5);
ind5 = find(mag>=5);
% mark into 8 classes


%% plot pmtm PSD
% fig 1, class 1~3
f1.fig=figure;
f1.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 3;     % number of stations
ncol = 3;   % number of components

for isub = 1: nrow*ncol
    f1.ax(isub) = subplot(nrow,ncol,isub);
end

if ~isempty(ind1)
    %%% class 1 N component
    subplot(nrow,ncol,1)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind1)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag<-1'),'FontSize',12,'unit','normalized');
    
    %%% class 1 E component
    subplot(nrow,ncol,2)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind1)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag<-1'),'FontSize',12,'unit','normalized');
    
    %%% class 1 U component
    subplot(nrow,ncol,3)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind1)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag<-1'),'FontSize',12,'unit','normalized');
end

if ~isempty(ind10)
    %%% class 2 N component
    subplot(nrow,ncol,4)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind10)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=-1&<0'),'FontSize',12,'unit','normalized');
    
    %%% class 2 E component
    subplot(nrow,ncol,5)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind10)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=-1&<0'),'FontSize',12,'unit','normalized');
    
    %%% class 2 U component
    subplot(nrow,ncol,6)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind10)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=-1&<0'),'FontSize',12,'unit','normalized');
end

if ~isempty(ind01)
    %%% class 3 N component
    subplot(nrow,ncol,7)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind01)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=0&<1'),'FontSize',12,'unit','normalized');
    
    %%% class 3 E component
    subplot(nrow,ncol,8)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind01)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=0&<1'),'FontSize',12,'unit','normalized');
    
    %%% class 3 U component
    subplot(nrow,ncol,9)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind01)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=0&<1'),'FontSize',12,'unit','normalized');
end


% fig 2, class 4~6
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 3;     % number of stations
ncol = 3;   % number of components

for isub = 1: nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

if ~isempty(ind12)
    %%% class 4 N component
    subplot(nrow,ncol,1)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind12)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=1&<2'),'FontSize',12,'unit','normalized');
    
    %%% class 4 E component
    subplot(nrow,ncol,2)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind12)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=1&<2'),'FontSize',12,'unit','normalized');
    
    %%% class 4 U component
    subplot(nrow,ncol,3)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind12)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=1&<2'),'FontSize',12,'unit','normalized');
end

% keyboard
if ~isempty(ind23)
    %%% class 5 N component
    subplot(nrow,ncol,4)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind23)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=2&<3'),'FontSize',12,'unit','normalized');
    
    %%% class 5 E component
    subplot(nrow,ncol,5)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind23)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=2&<3'),'FontSize',12,'unit','normalized');
    
    %%% class 5 U component
    subplot(nrow,ncol,6)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind23)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=2&<3'),'FontSize',12,'unit','normalized');
end

if ~isempty(ind34)
    %%% class 6 N component
    subplot(nrow,ncol,7)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind34)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=3&<4'),'FontSize',12,'unit','normalized');
    
    %%% class 6 E component
    subplot(nrow,ncol,8)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind34)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=3&<4'),'FontSize',12,'unit','normalized');
    
    %%% class 6 U component
    subplot(nrow,ncol,9)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind34)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=3&<4'),'FontSize',12,'unit','normalized');
end


% fig 3, class 7~8
f3.fig=figure;
f3.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 11;   % maximum height allowed is 11 inches
set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 2;     % number of stations
ncol = 3;   % number of components

for isub = 1: nrow*ncol
    f1.ax(isub) = subplot(nrow,ncol,isub);
end

if ~isempty(ind45)
    %%% class 7 N component
    subplot(nrow,ncol,1)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind45)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=4&<5'),'FontSize',12,'unit','normalized');
    
    %%% class 7 E component
    subplot(nrow,ncol,2)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind45)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=4&<5'),'FontSize',12,'unit','normalized');
    
    %%% class 7 U component
    subplot(nrow,ncol,3)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind45)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=4&<5'),'FontSize',12,'unit','normalized');
end

if ~isempty(ind5)
    %%% class 8 N component
    subplot(nrow,ncol,4)
    plt_psd_P_S_N_func(pmft,pmeq1NpsA,pmeq1NssA,pmeq1NnoiA,Fs,ind5)
    ax = gca;
    text(ax,0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=5'),'FontSize',12,'unit','normalized');
    
    %%% class 8 E component
    subplot(nrow,ncol,5)
    plt_psd_P_S_N_func(pmft,pmeq1EpsA,pmeq1EssA,pmeq1EnoiA,Fs,ind5)
    ax = gca;
    lgd = legend(ax,{'Noi bf ot','P sig','S sig','Noi med.','P med.','S med.'},'location','best',...
        'fontsize',10);
    olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
    set(olgd, 'Markersize', 12);
    text(ax,0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=5'),'FontSize',12,'unit','normalized');
    
    %%% class 8 U component
    subplot(nrow,ncol,6)
    plt_psd_P_S_N_func(pmft,pmeq1UpsA,pmeq1UssA,pmeq1UnoiA,Fs,ind5)
    ax = gca;
    text(ax,0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
    text(ax,0.1,0.1,strcat('mag>=5'),'FontSize',12,'unit','normalized');
end

% keyboard

% f2.fig=figure;
% f2.fig.Renderer='Painters';
% widin = 12;  % maximum width allowed is 8.5 inches
% htin = 4;   % maximum height allowed is 11 inches
% set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% 
% nrow = 1;     % number of stations
% ncol = 3;   % number of components
% 
% for isub = 1: nrow*ncol
%     f2.ax(isub) = subplot(nrow,ncol,isub);
% end
% 
% %%% N component
% subplot(nrow,ncol,1)
% grid on
% box on
% % noise before origin
% semilogx(pmft,10*log10(pmeq1Nnoi(:,1)),'-','color',[0.7 0.7 0.7],'linew',1); hold on
% % P signal
% semilogx(pmft,10*log10(pmeq1Nps(:,1)),'-','color',[0.3 0.3 0.3],'linew',1);
% % S signal
% semilogx(pmft,10*log10(pmeq1Nss(:,1)),'-','color','k','linew',1);
% 
% % tremor
% %         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% % ambient noise
% %         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);
% 
% xticks([1e-1, 1e0, 1e1]);
% xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB/Hz)');
% 
% text(0.1,0.95,strcat('N'),'FontSize',12,'unit','normalized');
% %          lgd = legend(ax,{'Noi bf ot','P sig','S sig'},'location','best',...
% %                       'fontsize',10);
% %          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
% %          set(olgd, 'Markersize', 12);
% hold off;
% 
% 
% %%% E component
% subplot(nrow,ncol,2)
% grid on
% box on
% % noise before origin
% semilogx(pmft,10*log10(pmeq1Enoi(:,1)),'-','color',[0.7 0.7 0.7],'linew',1); hold on
% % P signal
% semilogx(pmft,10*log10(pmeq1Eps(:,1)),'-','color',[0.3 0.3 0.3],'linew',1);
% % S signal
% semilogx(pmft,10*log10(pmeq1Ess(:,1)),'-','color','k','linew',1);
% 
% % tremor
% %         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% % ambient noise
% %         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);
% 
% xticks([1e-1, 1e0, 1e1]);
% xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
% xlabel('Frequency (Hz)');
% 
% text(0.1,0.95,strcat('E'),'FontSize',12,'unit','normalized');
% lgd = legend({'Noi bf ot','P sig','S sig'},'location','best',...
%     'fontsize',10);
% olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
% set(olgd, 'Markersize', 12);
% hold off;
% 
% 
% %%% U component
% subplot(nrow,ncol,3)
% grid on
% box on
% % noise before origin
% semilogx(pmft,10*log10(pmeq1Unoi(:,1)),'-','color',[0.7 0.7 0.7],'linew',1); hold on
% % P signal
% semilogx(pmft,10*log10(pmeq1Ups(:,1)),'-','color',[0.3 0.3 0.3],'linew',1);
% % S signal
% semilogx(pmft,10*log10(pmeq1Uss(:,1)),'-','color','k','linew',1);
% 
% % tremor
% %         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% % ambient noise
% %         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);
% 
% xticks([1e-1, 1e0, 1e1]);
% xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
% xlabel('Frequency (Hz)');
% 
% text(0.1,0.95,strcat('U'),'FontSize',12,'unit','normalized');
% %          lgd = legend({'Noi bf ot','P sig','S sig'},'location','best',...
% %                       'fontsize',10);
% %          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
% %          set(olgd, 'Markersize', 12);
% hold off;




















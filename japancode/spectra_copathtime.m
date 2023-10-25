function [pmeq1Ess,pmeq1Nss,pmeq1Uss,pmeq1Eps,pmeq1Nps,pmeq1Ups,pmft,f1,f2]=spectra_copathtime(...
    eq1Eps,eq1Nps,eq1Ups,eq1Ess,eq1Nss,eq1Uss,eq1Enoi,eq1Nnoi,eq1Unoi,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is to read in the P and S signal and noise window data of E
% N and U component data of hinet for each single target event that has co-path
% and-time tremor 
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/24
% Last modified date:   2020/07/25
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
% nfft = pow2(nextpow2(size(eq1Ess,1))-1);
nfft = pow2(nextpow2(size(eq1Ess,1)));
window = hann(size(eq1Ess,1));
Fs = sps;
[pdeq1Ess,pdft] = periodogram(eq1Ess,window,nfft,Fs);  % psd is done to each col, i.e., each station
[pdeq1Nss,~] = periodogram(eq1Nss,window,nfft,Fs);
[pdeq1Uss,~] = periodogram(eq1Uss,window,nfft,Fs);

[pdeq1Eps,~] = periodogram(eq1Eps,window,nfft,Fs);  % psd is done to each col, i.e., each station
[pdeq1Nps,~] = periodogram(eq1Nps,window,nfft,Fs);
[pdeq1Ups,~] = periodogram(eq1Ups,window,nfft,Fs);

% pmtm, P and S signal
nw = 4;
[pmeq1Ess,pmft] = pmtm(eq1Ess,nw,nfft,Fs);
[pmeq1Nss,~] = pmtm(eq1Nss,nw,nfft,Fs);
[pmeq1Uss,~] = pmtm(eq1Uss,nw,nfft,Fs);

[pmeq1Eps,~] = pmtm(eq1Eps,nw,nfft,Fs);
[pmeq1Nps,~] = pmtm(eq1Nps,nw,nfft,Fs);
[pmeq1Ups,~] = pmtm(eq1Ups,nw,nfft,Fs);


%%% EQ 1, noise around origin
% periodogram
% nfft = pow2(nextpow2(size(eq1Enoi,1))-1);
nfft = pow2(nextpow2(size(eq1Enoi,1)));
window = hann(size(eq1Enoi,1));
Fs = sps;
[pdeq1Enoi,~] = periodogram(eq1Enoi,window,nfft,Fs);
[pdeq1Nnoi,~] = periodogram(eq1Nnoi,window,nfft,Fs);
[pdeq1Unoi,~] = periodogram(eq1Unoi,window,nfft,Fs);

% pmtm
nw = 4;
[pmeq1Enoi,~] = pmtm(eq1Enoi,nw,nfft,Fs);
[pmeq1Nnoi,~] = pmtm(eq1Nnoi,nw,nfft,Fs);
[pmeq1Unoi,~] = pmtm(eq1Unoi,nw,nfft,Fs);


%% plot periodgram PSD
f1.fig=figure;
f1.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 1;     % number of stations
ncol = 3;   % number of components

for isub = 1: nrow*ncol
    f1.ax(isub) = subplot(nrow,ncol,isub);
end

%%% N component
subplot(nrow,ncol,1)
box on
% % noise before origin
% semilogx(pdft,10*log10(pdeq1Nnoi(:,1)),'-','color',[0.7 0.7 0.7],'linew',1); hold on; grid on;
% % P signal
% semilogx(pdft,10*log10(pdeq1Nps(:,1)),'-','color',[0.3 0.3 0.3],'linew',1);
% % S signal
% semilogx(pdft,10*log10(pdeq1Nss(:,1)),'-','color','k','linew',1);

% noise before origin
semilogx(pdft,10*log10(pdeq1Nnoi(:,1)),'-','color',[0.4 0.4 0.4],'linew',1.5); hold on; grid on;
% P signal
semilogx(pdft,10*log10(pdeq1Nps(:,1)),'-','color','r','linew',1.5);
% S signal
semilogx(pdft,10*log10(pdeq1Nss(:,1)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pdft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pdft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
ylim([-320,-120]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');

text(0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
%          lgd = legend(ax,{'Noi bf ot','P sig','S sig'},'location','best',...
%                       'fontsize',10);
%          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
%          set(olgd, 'Markersize', 12);
hold off;


%%% E component
subplot(nrow,ncol,2)
box on
% noise before origin
semilogx(pdft,10*log10(pdeq1Enoi(:,1)),'-','color',[0.4 0.4 0.4],'linew',1.5); hold on; grid on;
% P signal
semilogx(pdft,10*log10(pdeq1Eps(:,1)),'-','color','r','linew',1.5);
% S signal
semilogx(pdft,10*log10(pdeq1Ess(:,1)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pdft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pdft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
ylim([-320,-120]);
xlabel('Frequency (Hz)');

text(0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
lgd = legend({'Noi bf ot','P sig','S sig'},'location','best',...
    'fontsize',10);
olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
set(olgd, 'Markersize', 12);
hold off;


%%% U component
subplot(nrow,ncol,3)
box on
% noise before origin
semilogx(pdft,10*log10(pdeq1Unoi(:,1)),'-','color',[0.4 0.4 0.4],'linew',1.5); hold on; grid on;
% P signal
semilogx(pdft,10*log10(pdeq1Ups(:,1)),'-','color','r','linew',1.5);
% S signal
semilogx(pdft,10*log10(pdeq1Uss(:,1)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pdft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pdft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
ylim([-320,-120]);
xlabel('Frequency (Hz)');

text(0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
%          lgd = legend({'Noi bf ot','P sig','S sig'},'location','best',...
%                       'fontsize',10);
%          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
%          set(olgd, 'Markersize', 12);
hold off;


%% plot pmtm PSD
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 1;     % number of stations
ncol = 3;   % number of components

for isub = 1: nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

%%% N component
subplot(nrow,ncol,1)
box on
% noise before origin
semilogx(pmft,10*log10(pmeq1Nnoi(:,1)),'-','color',[0.4 0.4 0.4],'linew',1.5); hold on; grid on;
% P signal
semilogx(pmft,10*log10(pmeq1Nps(:,1)),'-','color','r','linew',1.5);
% S signa
semilogx(pmft,10*log10(pmeq1Nss(:,1)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');

text(0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
%          lgd = legend(ax,{'Noi bf ot','P sig','S sig'},'location','best',...
%                       'fontsize',10);
%          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
%          set(olgd, 'Markersize', 12);
hold off;


%%% E component
subplot(nrow,ncol,2)
box on
% noise before origin
semilogx(pmft,10*log10(pmeq1Enoi(:,1)),'-','color',[0.4 0.4 0.4],'linew',1.5); hold on; grid on;
% P signal
semilogx(pmft,10*log10(pmeq1Eps(:,1)),'-','color','r','linew',1.5);
% S signa
semilogx(pmft,10*log10(pmeq1Ess(:,1)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
xlabel('Frequency (Hz)');

text(0.8,0.95,strcat('E'),'FontSize',12,'unit','normalized');
lgd = legend({'Noi bf ot','P sig','S sig'},'location','best',...
    'fontsize',10);
olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
set(olgd, 'Markersize', 12);
hold off;


%%% U component
subplot(nrow,ncol,3)
box on
% noise before origin
semilogx(pmft,10*log10(pmeq1Unoi(:,1)),'-','color',[0.4 0.4 0.4],'linew',1.5); hold on; grid on;
% P signal
semilogx(pmft,10*log10(pmeq1Ups(:,1)),'-','color','r','linew',1.5);
% S signa
semilogx(pmft,10*log10(pmeq1Uss(:,1)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
xlabel('Frequency (Hz)');

text(0.8,0.95,strcat('U'),'FontSize',12,'unit','normalized');
%          lgd = legend({'Noi bf ot','P sig','S sig'},'location','best',...
%                       'fontsize',10);
%          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
%          set(olgd, 'Markersize', 12);
hold off;




















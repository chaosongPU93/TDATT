function plt_psd_P_S_N_func(pmft,pmeqps,pmeqss,pmeqnoi,Fs,ind)

box on
% noise before origin, light gray
%semilogx(pmft,10*log10(pmeqnoi(:,ind)),'-','color',[0.7 0.7 0.7],'linew',1); hold on
semilogx(pmft,10*log10(pmeqnoi(:,ind)),'-','color',[0.8 0.8 0.8],'linew',0.5); hold on; grid on;
% P signal, light red
%semilogx(pmft,10*log10(pmeqps(:,ind)),'-','color',[0.3 0.3 0.3],'linew',1);
semilogx(pmft,10*log10(pmeqps(:,ind)),'-','color',[0.95 0.5 0.5],'linew',0.5);
% S signal, light blue
%semilogx(pmft,10*log10(pmeqss(:,ind)),'-','color','k','linew',1);
semilogx(pmft,10*log10(pmeqss(:,ind)),'-','color',[0 1 1],'linew',0.5);

% next plot their median
semilogx(pmft,10*log10(median(pmeqnoi(:,ind),2)),'-','color',[0.4,0.4,0.4],'linew',1.5);
semilogx(pmft,10*log10(median(pmeqps(:,ind),2)),'-','color','r','linew',1.5);
semilogx(pmft,10*log10(median(pmeqss(:,ind),2)),'-','color','b','linew',1.5);

% tremor
%         semilogx(pmft,10*log10(pdtmrN(:,i)),'--','color','k','linew',1);
% ambient noise
%         semilogx(pmft,10*log10(pdnoiN(:,i)),':','color','k','linew',1);

xticks([1e-1, 1e0, 1e1]);
xlim([1e-1, Fs/2]);
% ylim([-320,-120]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');

% text(0.8,0.95,strcat('N'),'FontSize',12,'unit','normalized');
%          lgd = legend(ax,{'Noi bf ot','P sig','S sig'},'location','best',...
%                       'fontsize',10);
%          olgd = findobj(lgd, 'type', 'line'); %// objects of legend of type line
%          set(olgd, 'Markersize', 12);
hold off;

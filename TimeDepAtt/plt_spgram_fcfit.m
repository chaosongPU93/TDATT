function [f,fitobj,gof,outp,fgmaxfit] = plt_spgram_fcfit(tets,fgmax,fran,yran)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f,fitobj,gof,outp,fgmaxfit] = plt_spgram_fcfit(tets,fgmax,fran,yran)
% This function is plot the automatically-picked corner frequencies (or should
% be called peaks in the spectrogram within a frequency range). A range of freq
% that has define some percent of the max value is viewed as an error bar. 
% A robust linear fit using the bisquare scheme is also plotted. Basically uses
% the outputs from 'plt_spgram_of_bursts_norm_zoom'
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/09
% Last modified date:   2021/11/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
f.fig=figure;
[scrsz, res] = pixelperinch(1);

ax = gca; 
box on; 
grid on; 
hold on
errorbar(ax,tets,fgmax,fgmax-fran(:,1),fran(:,2)-fgmax,'vertical','o','markersize',4.5,'color',...
         [.7 .7 .7],'linewidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',5);
% linear robust least square
[fitobj,gof,outp,fgmaxfit] = linear_bisquare_fit_free(tets,fgmax);
plot(ax,tets,fgmaxfit,'-','linewidth',3,'color','r');
coef = coeffvalues(fitobj);
slope = coef(1);
disp(slope);
disp(fgmaxfit(end)-fgmaxfit(1));
disp(fgmaxfit(1));
disp((fgmaxfit(end)-fgmaxfit(1))/fgmaxfit(1));
text(ax,0.6, 0.95, sprintf('Slope: %.1e',slope),'fontsize',10,'unit','normalized',...
        'horizontalalignment','left');
fcchg = fgmaxfit(end)-fgmaxfit(1);
text(ax,0.6, 0.85, sprintf('Change: %.2f/%.2f=%.1f%%',fcchg,fgmaxfit(1),...
  fcchg/fgmaxfit(1)*100),'fontsize',10,'unit','normalized',...
        'horizontalalignment','left');
ax.YScale = 'log';
xlim(ax,[min(tets) max(tets)]);
if ~isempty(yran)
  ylim(ax,yran);
end
ylabel(ax, 'Frequency (Hz)','fontsize',12);
xlabel(ax,strcat({'Time (s)'}),'fontsize',12);
axes('xlim',ax.XLim,'ylim',ax.YLim,'color','none','YAxisLocation','right','XTick',[],...
  'YScale','log','YTick',ax.YTick);


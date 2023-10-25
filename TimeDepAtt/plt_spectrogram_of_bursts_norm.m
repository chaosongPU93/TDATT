function [f] = plt_spectrogram_of_bursts_norm(years,stas,iets,spgft,spgall,tall,minfnorm,maxfnorm,yran,cran,...
                                              filtflag,filtsigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to plot the amplitude spectra of all the tremor burst
% windows of each ETS, ignoring information of time. Each spectrum is
% normalized with the mean amplitude within a frequency range. Each trace
% is shaded in gray, while the 'median' could be reinforced by the
% strengthening of color due to the stacking
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/08
% Last modified date:   2021/11/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f.fig = figure;
%   f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 10.6;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nsta = size(stas,1);
nrow = nsta;
ncol = 1;

for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
  f.ax(isub).Box = 'on';
  %     grid(f.ax(isub),'on');
end

%%% reposition
for ista = 1: nsta
  set(f.ax(ista), 'position', [ 0.08, 0.08+(0.12/(nsta-1)+0.75/nsta)*(nsta-ista), 0.78, 0.8/nsta]);
end
%     set(f.ax(1), 'position', [ 0.08, 0.08, 0.8, 0.8]);
%   set(f.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
%   set(f.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);

isub = 0;
pltets = spgall{iets,1};
tets = tall{iets,1};
for ista = 1: nsta
  
  isub = isub+1;
  ax = f.ax(isub);
  hold(ax,'on');
  
  pltsta = squeeze(pltets(:,:,ista));
  pltsta = sqrt(pltsta);  % convert power to amplitude
  
  % normalize the spectrogram according to the mean amplitude within a freq range, note that each
  % col of the entire array is the FFT of each overlapping window, where there is no averaging
  %     minfreq = 2; maxfreq = 4;
  %     minfreq = 1.5; maxfreq = 5;
%   minfnorm = 2.5; maxfnorm = 4;
  [~,indmin] = min(abs(spgft-minfnorm));
  [~,indmax] = min(abs(spgft-maxfnorm));
  normlzer = mean(pltsta(indmin:indmax,:),1);
  pltstanorm=pltsta;
  %     keyboard
  for k = 1: size(pltsta, 2)  %here k is number of overlapping subwindow for FFT, larger than number of bursts
    pltstanorm(:,k)=pltstanorm(:,k)/normlzer(k);
  end   % all overlapping subwindows for FFT
  
  if filtflag == 1
    %%% smooth with a gaussian filter, to make it more continuous in time, requires more smoothing
    %%% in time, ie, col, dimension
    %     pltstanorm = imgaussfilt(pltstanorm,[5*size(pltstanorm,1)/size(pltstanorm,2) 5]);
    pltstanorm = imgaussfilt(pltstanorm,filtsigma);
  end
  
  %%% freq range to plot
  minfplt = yran(1); maxfplt = yran(2); %minfreq=2.;
  [~,indmin] = min(abs(spgft-minfplt));
  [~,indmax] = min(abs(spgft-maxfplt));
  surf(ax,tets,spgft(indmin:indmax),log10(pltstanorm(indmin:indmax, :)),'EdgeColor','none');

%   surf(ax,tets,spgft,log10(pltstanorm),'EdgeColor','none');
%     surf(ax,tets,spgft,pltstanorm,'EdgeColor','none');
  
  %     ylim(ax,[0.1 20]);
  %     ylim(ax,[1 10]);
  %     yticks(ax,[0.1 1 10]);
  ylim(ax,yran);
  yticks(ax,1:1:10);
  xlim(ax,[min(tets) max(tets)]);
  ax.YScale = 'log';
  %     ax.YDir = 'reverse';
  ylabel(ax, 'Frequency (Hz)','fontsize',12);
  if ista == nsta
    xlabel(ax,strcat({'Time (s)'}),'fontsize',12);
  end
  %     colormap(ax,'jet');
  %     colormap(ax,'redblue');
  oldc = colormap(ax,'kelicol');
  newc = flipud(oldc);
  colormap(ax,newc);
  c=colorbar(ax,'East');
  pos = ax.Position;
  c.Position = [pos(1)+pos(3)+0.03, pos(2), 0.02, pos(4)];
  c.Label.String = 'log_{10}(Amplitude)/Hz';
  c.FontSize = 10;
  %     c.Label.String = 'Amplitude/Hz';
  if ~isempty(cran)
    caxis(ax,cran);
  end
%   cran = min(abs(caxis(ax)));
  %     caxis(ax,[-cran*0.5,cran*0.5]);
  %     caxis(ax,[-cran*0.8,cran*0.8]);
  view(ax,0,90);
  ax.TickDir = 'out';
  
  %     ax.GridLineStyle = '-';
  text(ax,0.04,0.9,num2str(years(iets)),'FontSize',10,'unit','normalized',...
    'horizontalalignment','left','backgroundColor','w','Margin',1,'color','k');
  text(ax,0.96,0.9,stas(ista, :),'FontSize',10,'unit','normalized',...
    'horizontalalignment','right','backgroundColor','w','Margin',1,'color','k');
  
  axes('xlim',ax.XLim,'ylim',ax.YLim,'color','none','YAxisLocation','right','XTick',[],...
    'YScale','log','YTick',1:1:10,'position', ax.Position);
  
end   % all stations

function [f] = plt_spectra_of_bursts(years,stas,pcft,pcall,xran,yran)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to plot the amplitude spectra of all the tremor burst 
% windows of each ETS, ignoring information of time. Each trace
% is shaded in gray, while the 'median' could be reinforced by the 
% strengthening of color due to the stacking.
% There is NO normalization, and the counterpart 'plt_spectra_of_bursts_norm'
% normalizes each spectrum with the mean amplitude within a frequency range.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/06/15
% Last modified date:   2022/06/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nets = length(years);
nsta = size(stas,1);
nrow = nets;
ncol = nsta;

for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
  f.ax(isub).Box = 'on';
  %     grid(f.ax(isub),'on');
end

pltxran = [0.08 0.98]; pltyran = [0.08 0.98];
pltxsep = 0.05; pltysep = 0.05; 
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

isub = 0;
for iets = 1: nets
  pltets = pcall{iets,1};

  for ista = 1: nsta
    
    isub = isub+1;
    ax = f.ax(isub);
%     hold(ax,'on');
    
    pltsta = squeeze(pltets(:,ista,:));
    pltsta = sqrt(pltsta);  % convert power to amplitude
    %note that for each single spectrum k (one col), it is already a 'mean' of several
    %overlapping windows
    for k = 1: size(pltsta, 2)  %here k is number of the burst window
      loglog(ax,pcft,pltsta(:,k),'color',[0 0 0  0.15])
      hold(ax,'on');
%       grid(ax,'on');  
    end   % all bursts

    pltstamed = median(pltsta,2);  % median of all cols, ie, bursts
    pltstaave = mean(pltsta,2);  % mean of all cols, ie, bursts
    loglog(ax,pcft,pltstaave,'r','linew',1);
    loglog(ax,pcft,pltstamed,'b','linew',1);
    
    loglog(ax,[1 10],[1 0.1],'g','linew',1.5);  % a reference line with a slope of -1
    
%     axis(ax, 'equal');
    xlim(ax,xran);
    ylim(ax,yran);
    xticks(ax,[0.1 1 10]);
    
    plot(ax,[1 1],ax.YLim,':','color',[0.8 0.8 0.8]);
    plot(ax,[2 2],ax.YLim,':','color',[0.8 0.8 0.8]);
    plot(ax,[4 4],ax.YLim,':','color',[0.8 0.8 0.8]);
    plot(ax,[8 8],ax.YLim,':','color',[0.8 0.8 0.8]);

%     ax.GridLineStyle = '-';
    text(ax,0.04,0.95,num2str(years(iets)),'FontSize',8,'unit','normalized',...
      'horizontalalignment','left');  %,'EdgeColor','k','Margin',1
    text(ax,0.96,0.95,stas(ista, :),'FontSize',8,'unit','normalized',...
      'horizontalalignment','right');

  end   % all stations
end   % all ets

ax = f.ax((nrow-1)*ncol+1);
hold(ax,'on');
xlabel(ax,'Frequency (Hz)','FontSize',12);
ylabel(ax,'Amplitude/Hz','FontSize',12);


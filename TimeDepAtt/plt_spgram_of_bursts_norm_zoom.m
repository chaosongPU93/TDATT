function [f,tets,fgmax,fgmaxfit,fran] = plt_spgram_of_bursts_norm_zoom(years,stas,iets,ista,...
          spgft,spgall,tall,minfnorm,maxfnorm,yran,cran,filtflag,filtsigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is similar to 'plt_spectrogram_of_bursts_norm', but to zoom
% at one station to get a better vertical/horizontal figure ratio
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

nrow = 1;
ncol = 1;

for isub = 1:nrow*ncol
  f.ax(isub) = subplot(nrow,ncol,isub);
  f.ax(isub).Box = 'on';
  %     grid(f.ax(isub),'on');
end

%%% reposition
set(f.ax(1), 'position', [ 0.08, 0.08, 0.8, 0.8]);

isub = 0;
pltets = spgall{iets,1};
tets = tall{iets,1};

%   ista = 4;
isub = isub+1;
ax = f.ax(isub);
hold(ax,'on');

pltsta = squeeze(pltets(:,:,ista));
pltsta = sqrt(pltsta);  % convert power to amplitude

% normalize the spectrogram according to the mean amplitude within a freq range, note that each
% col of the entire array is the FFT of each overlapping window, where there is no averaging
%     minfreq = 2; maxfreq = 4; %minfreq=2.;
% minfnorm = 1.5; maxfnorm = 5;
%     minfreq = 3; maxfreq = 4;
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
  
%   h = fspecial('gaussian',[],2);
%   pltstanorm = imfilter(pltstanorm,h);
end

%%% freq range to plot
minfplt = yran(1); maxfplt = yran(2); %minfreq=2.;
[~,indmin] = min(abs(spgft-minfplt));
[~,indmax] = min(abs(spgft-maxfplt));
surf(ax,tets,spgft(indmin:indmax),log10(pltstanorm(indmin:indmax, :)),'EdgeColor','none');

%     imagesc(ax, [tets(1) tets(end)], [spgft(indmin), spgft(indmax)],log10(pltstanorm(indmin:indmax, :)));

% surf(ax,tets,spgft,log10(pltstanorm),'EdgeColor','none');
%     surf(ax,tets,spgft,pltstanorm,'EdgeColor','none');

ylim(ax,yran);
yticks(ax,1:1:10);
xlim(ax,[min(tets) max(tets)]);
ax.YScale = 'log';
%     ax.YDir = 'reverse';
ylabel(ax, 'Frequency (Hz)','fontsize',10);
xlabel(ax,strcat({'Time (s)'}),'fontsize',10);
%     colormap(ax,'jet');
%     colormap(ax,'redblue');
%     colormap(ax,'maroonnavy');
oldc = colormap(ax,'kelicol');
newc = flipud(oldc);
colormap(ax,newc);
c=colorbar(ax,'East');
pos = ax.Position;
c.Position = [pos(1)+pos(3)+0.04, pos(2), 0.02, pos(4)];
c.Label.String = 'log_{10}(Amplitude)/Hz';
%     c.Label.String = 'Amplitude/Hz';
%     if iets == 1
%       caxis(ax,[-1,2.5]);
%     elseif iets == 2
%       caxis(ax,[-1,1.7]);
%     elseif iets == 3
%       caxis(ax,[-1,1.7]);
%     end
%     caxis(ax,[-0.7,0.7]);
% cran = min(abs(caxis(ax)));
%     caxis(ax,[-cran*0.5,cran*0.5]);
%     caxis(ax,[-cran*0.8,cran*0.8]);
caxis(ax,cran);

ax.TickDir = 'out';
%     view(ax,0,90);
view(2);

%     ax.GridLineStyle = '-';
text(ax,0.04,0.9,num2str(years(iets)),'FontSize',10,'unit','normalized',...
  'horizontalalignment','left','backgroundColor','w','Margin',1,'color','k');
text(ax,0.96,0.9,stas(ista, :),'FontSize',10,'unit','normalized',...
  'horizontalalignment','right','backgroundColor','w','Margin',1,'color','k');

axes('xlim',ax.XLim,'ylim',ax.YLim,'color','none','YAxisLocation','right','XTick',[],...
  'YScale','log','YTick',1:1:10,'position', [ 0.08, 0.08, 0.8, 0.8]);


%%%if Gaussian filtering is used, then it is worthwhile to find the peaks, range, and do regression
if filtflag
  
%%% find the maxima of the each col. between the freq of interest
if ista == 3
  minfpeak = 3.3; maxfpeak = 4.27; %minfreq=2.;
elseif ista == 4
  minfpeak = 2.9; maxfpeak = 3.85; %minfreq=2.;
elseif ista == 2
  minfpeak = 2.7; maxfpeak = 3.25; %minfreq=2.;
end
plot3(ax,ax.XLim, [minfpeak minfpeak],10*ones(2,1),'w--','linewidth',2);
plot3(ax,ax.XLim, [maxfpeak maxfpeak],10*ones(2,1),'w--','linewidth',2);
[~,indmin] = min(abs(spgft-minfpeak));
[~,indmax] = min(abs(spgft-maxfpeak));
[valmax, fgmaxind] = max(pltstanorm(indmin:indmax,:), [], 1);
fgmaxind = fgmaxind+indmin-1; % convert to global index
fgmax = spgft(fgmaxind);
%     plot(ax,tets, fmax,'o-','linewidth',1,'color','k','markers',1.5);
scatter3(ax, tets, fgmax, 10*ones(length(tets),1),30,'ko','filled');

% find the freq that has some percent of the global max value that is useful to define to range
fran = zeros(size(pltstanorm,2), 2);
buff = 20;  % allow some buffer in freq
for i = 1: size(pltstanorm, 2)
  
  [lpkval,lpkind] = findpeaks(pltstanorm(1:fgmaxind(i)-1,i),...
                                                      'MinPeakDistance', 4);
  lonlpkind = lpkind(end);
  buff = ceil((fgmaxind(i)-lonlpkind)*2/3);
  [~, franindlo] = min(abs(pltstanorm(fgmaxind(i)-buff:fgmaxind(i)-1,i)-0.95*valmax(i)));
  if isempty(franindlo)
    disp(i)
    break
  end
  franindlo = franindlo+fgmaxind(i)-buff-1; % convert to global index
  fran(i,1) = spgft(franindlo);
  
  [lpkval,lpkind] = findpeaks(pltstanorm(fgmaxind(i)+1: end,i),...
                                                      'MinPeakDistance', 4);
  hinlpkind = lpkind(1)+fgmaxind(i);
  buff = ceil((hinlpkind-fgmaxind(i))*2/3);
  [~, franindlo] = min(abs(pltstanorm(fgmaxind(i)+1:fgmaxind(i)+buff,i)-0.95*valmax(i)));
  if isempty(franindlo)
    disp(i)
    break
  end
  franindlo = franindlo+fgmaxind(i); % convert to global index
  fran(i,2) = spgft(franindlo);
                                                  
%   [~, franindlo] = min(abs(pltstanorm(fgmaxind(i)-20:fgmaxind(i),i)-0.95*valmax(i)));
%   if isempty(franindlo)
%     disp(i)
%     break
%   end
%   franindlo = franindlo+fgmaxind(i)-20-1; % convert to global index
%   fran(i,1) = spgft(franindlo);
%   [~, franindhi] = min(abs(pltstanorm(fgmaxind(i):fgmaxind(i)+buff,i)-0.95*valmax(i)));
%   if isempty(franindhi)
%     disp(i)
%     break
%   end
%   franindhi = franindhi+fgmaxind(i)-1; % convert to global index
%   fran(i,2) = spgft(franindhi);
end
       
% linear robust least square
[~,~,~,fgmaxfit] = linear_bisquare_fit_free(tets,fgmax);
% coef = coeffvalues(fitobj);
% slope = coef(1);
% intcpt = coef(2);
plot3(ax,tets,fgmaxfit,10*ones(length(fgmaxfit)),'-','linewidth',3,'color',[.6 .6 .6]);

if ista == 3
  minfpeak = 4.35; maxfpeak = 5; %minfreq=2.;
  plot3(ax,ax.XLim, [minfpeak minfpeak],10*ones(2,1),'k--','linewidth',2);
  plot3(ax,ax.XLim, [maxfpeak maxfpeak],10*ones(2,1),'k--','linewidth',2);
  [~,indmin] = min(abs(spgft-minfpeak));
  [~,indmax] = min(abs(spgft-maxfpeak));
  [valmax, fgmaxind] = max(pltstanorm(indmin:indmax,:), [], 1);
  fgmaxind = fgmaxind+indmin-1; % convert to global index
  fgmax(:,2) = spgft(fgmaxind);
  %     plot(ax,tets, fmax,'o-','linewidth',1,'color','k','markers',1.5);
  scatter3(ax, tets, fgmax(:,2), 10*ones(length(tets),1),30,'ks','filled');
  
  % linear robust least square
  [~,~,~,fgmaxfit(:,2)] = linear_bisquare_fit_free(tets,fgmax(:,2));
  plot3(ax,tets,fgmaxfit(:,2),10*ones(length(fgmaxfit(:,2))),'-','linewidth',3,'color',[.6 .6 .6]);

  % find the freq that has some percent of the global max value that is useful to define to range
%   fran = zeros(size(pltstanorm,2), 2);
  buff = 20;  % allow some buffer in freq
  for i = 1: size(pltstanorm, 2)
    [lpkval,lpkind] = findpeaks(pltstanorm(1:fgmaxind(i)-1,i),...
      'MinPeakDistance', 4);
    lonlpkind = lpkind(end);
    buff = ceil((fgmaxind(i)-lonlpkind)*2/3);
    [~, franindlo] = min(abs(pltstanorm(fgmaxind(i)-buff:fgmaxind(i)-1,i)-0.95*valmax(i)));
    if isempty(franindlo)
      disp(i)
      break
    end
    franindlo = franindlo+fgmaxind(i)-buff-1; % convert to global index
    fran(i,3) = spgft(franindlo);
    
    [lpkval,lpkind] = findpeaks(pltstanorm(fgmaxind(i)+1: end,i),...
      'MinPeakDistance', 4);
    hinlpkind = lpkind(1)+fgmaxind(i);
    buff = ceil((hinlpkind-fgmaxind(i))*2/3);
    [~, franindlo] = min(abs(pltstanorm(fgmaxind(i)+1:fgmaxind(i)+buff,i)-0.95*valmax(i)));
    if isempty(franindlo)
      disp(i)
      break
    end
    franindlo = franindlo+fgmaxind(i); % convert to global index
    fran(i,4) = spgft(franindlo);

%     [~, franindlo] = min(abs(pltstanorm(fgmaxind(i)-buff:fgmaxind(i),i)-0.95*valmax(i)));
%     if isempty(franindlo)
%       disp(i)
%       break
%     end
%     franindlo = franindlo+fgmaxind(i)-buff-1; % convert to global index
%     fran(i,3) = spgft(franindlo);
%     [~, franindhi] = min(abs(pltstanorm(fgmaxind(i):fgmaxind(i)+buff,i)-0.95*valmax(i)));
%     if isempty(franindhi)
%       disp(i)
%       break
%     end
%     franindhi = franindhi+fgmaxind(i)-1; % convert to global index
%     fran(i,4) = spgft(franindhi);
  end

end

else
  fgmax = []; fgmaxfit = []; fran = [];
end


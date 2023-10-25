function [f] = plt_burst_diffreg(ymax,burstcell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to plot the selected tremor bursts in time as gray bar
% in the same ETS episode for different regions
%
% Since the number of regions to be ploted is variable, so here the input
% 'burstcell' is a cell array with each element being the bursts of each
% region
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/10/26
% Last modified date:   2021/10/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[scrsz, res] = pixelperinch(1);

nrow = size(burstcell,1);
ncol = 1;

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 1.5*nrow;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);


for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
end

for i = 1: nrow
  ax = f.ax(i);
  hold(ax,'on');
  
  burst = burstcell{i};
  t = zeros(size(burst,1),2);
  if ~isempty(burst)
    if burst(1,1) < 2004*1000
      t(:,1) = (burst(:,1)-2003060)+burst(:,2)./(3600.*24);
      t(:,2) = (burst(:,1)-2003060)+burst(:,3)./(3600.*24);
      text(f.ax(1),0.92,0.9,'2003','FontSize',11,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
      xlabel(f.ax(nrow), 'Time (day) since Mar. 1, 2003 (2003060)');
    elseif burst(1,1) > 2004*1000 && burst(1,1)< 2005*1000
      t(:,1) = (burst(:,1)-2004194)+burst(:,2)./(3600.*24);
      t(:,2) = (burst(:,1)-2004194)+burst(:,3)./(3600.*24);
      text(f.ax(1),0.92,0.9,'2004','FontSize',11,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
      xlabel(f.ax(nrow), 'Time (day) since Jul. 12, 2004 (2004194)');
    elseif burst(1,1) > 2005*1000
      t(:,1) = (burst(:,1)-2005254)+burst(:,2)./(3600.*24);
      t(:,2) = (burst(:,1)-2005254)+burst(:,3)./(3600.*24);
      text(f.ax(1),0.92,0.9,'2005','FontSize',11,'unit','normalized','horizontalalignment','left',...
        'EdgeColor','k','Margin',2);
      xlabel(f.ax(nrow), 'Time (day) since Sep. 11, 2005 (2005254)');
    end
    
    for j = 1: size(t,1)
      patarea = [t(j,1) -3e+4;
        t(j,2) -3e+4;
        t(j,2) 3e+4;
        t(j,1) 3e+4;
        t(j,1) -3e+4];
      patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
  end
  
%   ylabel(ax, 'Cumulative number of detections');
  xlim(ax,f.ax(1).XLim);
  ylim(ax,[0, ymax]);
  % ytickformat(ax,'%e');
  hold(ax,'off');
  
end



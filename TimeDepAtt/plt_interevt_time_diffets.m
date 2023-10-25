function [f] = plt_interevt_time_diffets(hfinter03,hfinter04,hfinter05,ymax,ttol1,ttol2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_interevt_time_diffets(hfinter03,hfinter04,hfinter05,ymax,ttol1,ttol2)
% This function is to plot the 
% separation in time between itself and its preceding detection,
% i.e., inter-event time, for all 3 ETS episodes
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/26
% Last modified date:   2022/03/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

defval('ymax',2e-3);
defval('ttol1',[]);
defval('ttol2',[]);

[scrsz, res] = pixelperinch(1);

f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 3;
ncol = 1;

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
end

msizehf = 2;

ax = f.ax(1);
hold(ax,'on'); 
scatter(ax, hfinter03(:,1), hfinter03(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
if ~isempty(ttol1)
  plot(ax, ax.XLim, [ttol1(1) ttol1(1)], 'k--');
end
if ~isempty(ttol2)
  plot(ax, ax.XLim, [ttol2(1) ttol2(1)], 'b--');
end
text(ax,0.92,0.9,'2003','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
ylim(ax,[0, ymax]);
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); 
scatter(ax, hfinter04(:,1), hfinter04(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
if ~isempty(ttol1)
  plot(ax, ax.XLim, [ttol1(2) ttol1(2)], 'k--');
end
if ~isempty(ttol2)
  plot(ax, ax.XLim, [ttol2(2) ttol2(2)], 'b--');
end
text(ax,0.92,0.9,'2004','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since Jul. 12, 2004 (2004194)');
ylim(ax,[0, ymax]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); 
scatter(ax, hfinter05(:,1), hfinter05(:,2), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
if ~isempty(ttol1)
  plot(ax, ax.XLim, [ttol1(3) ttol1(3)], 'k--');
end
if ~isempty(ttol2)
  plot(ax, ax.XLim, [ttol2(3) ttol2(3)], 'b--');
end
text(ax,0.92,0.9,'2005','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Time (day) to the preceding detection');
xlabel(ax, 'Time (day) since Sep. 11, 2005 (2005254)');
ylim(ax,[0, ymax]);
hold(ax,'off');





function [f] = plt_burst_diffets(ymax,newt03,newt04,newt05,hfinter03,hfinter04,hfinter05)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f] = plt_burst_diffets(ymax,newt03,newt04,newt05,hfinter03,hfinter04,hfinter05)
% This function is to plot the selected tremor bursts in time as gray bar
% and the cumulative number of detections for the same region at different
% ETS episodes
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/26
% Last modified date:   2021/09/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

ax = f.ax(1);
hold(ax,'on');
if ~isempty(newt03)
    for j = 1: size(newt03,1)
        patarea = [newt03(j,1) -3e+4;
                   newt03(j,2) -3e+4;
                   newt03(j,2) 3e+4;
                   newt03(j,1) 3e+4;
                   newt03(j,1) -3e+4];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
end
stairs(ax, hfinter03(:,1), 1: length(hfinter03(:,1)),'b-');  %, 'MarkerEdgeColor', 'w')
text(ax,0.92,0.9,'2003','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% text(ax,0.06,0.8,strcat(num2str(round(perchf03)),'%'),'FontSize',12,'unit','normalized','horizontalalignment',...
%     'left');
ylabel(ax, 'Cumulative number of detections');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
ylim(ax,[0, ymax]);
% ytickformat(ax,'%e');
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); 
if ~isempty(newt04)
    for j = 1: size(newt04,1)
        patarea = [newt04(j,1) -3e+4;
                   newt04(j,2) -3e+4;
                   newt04(j,2) 3e+4;
                   newt04(j,1) 3e+4;
                   newt04(j,1) -3e+4];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
end
stairs(ax, hfinter04(:,1), (1: length(hfinter04(:,1)))','b-');  %, 'MarkerEdgeColor', 'w')
text(ax,0.92,0.9,'2004','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% text(ax,0.06,0.8,strcat(num2str(round(perchf04)),'%'),'FontSize',12,'unit','normalized','horizontalalignment',...
%     'left');
% ylabel(ax, 'Cumulative number of detections');
xlabel(ax, 'Time (day) since Jul. 12, 2004 (2004194)');
ylim(ax,[0, ymax]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); 
if ~isempty(newt05)
    for j = 1: size(newt05,1)
        patarea = [newt05(j,1) -3e+4;
                   newt05(j,2) -3e+4;
                   newt05(j,2) 3e+4;
                   newt05(j,1) 3e+4;
                   newt05(j,1) -3e+4];
        patch(ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.3,'edgecolor','none');
    end
end
stairs(ax, hfinter05(:,1), (1: length(hfinter05(:,1)))','b-');  %, 'MarkerEdgeColor', 'w')
text(ax,0.92,0.9,'2005','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
% text(ax,0.06,0.8,strcat(num2str(round(perchf05)),'%'),'FontSize',12,'unit','normalized','horizontalalignment',...
%     'left');
% ylabel(ax, 'Cumulative number of detections');
xlabel(ax, 'Time (day) since Sep. 11, 2005 (2005254)');
ylim(ax,[0, ymax]);
hold(ax,'off');




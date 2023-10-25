%   collect all the unused parts but potentially useful in the future of the codes in 'proj2021'

% Chao, 2021/05/27






%% from 'identify_pgc.m' and/or 'identify.m'
allrst_new = sortrows(allrst_new, [1, 8]);      % use the time of the strongest arrival

il03 = find(allrst_new(:,1) < 2004*1000);
allrst_new(il03,20) = (allrst_new(il03,1)-2003060)+allrst_new(il03,8)./(3600.*24);
il04 = find(allrst_new(:,1) < 2005*1000 & allrst_new(:,1) > 2004*1000);
allrst_new(il04,20) = (allrst_new(il04,1)-2004194)+allrst_new(il04,8)./(3600.*24);
il05 = find(allrst_new(:,1) > 2005*1000);
allrst_new(il05,20) = (allrst_new(il05,1)-2005254)+allrst_new(il05,8)./(3600.*24);

[allrst_sort, indori] = sortrows(allrst_new, [11, 12],'descend');      % use the time of the strongest arrival

%%%
%%% plot the energy ratio
f.fig = figure;
f.fig.Renderer = 'painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 8;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 4;
ncol = 1;

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
end

msizehf = 2;

ax = f.ax(1);
hold(ax,'on'); 
s1 = scatter(ax, allrst_new(:,end), allrst_new(:,11), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
s2 = scatter(ax, allrst_new(:,end), allrst_new(:,12), msizehf, 'filled','bo');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
% plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'2.5*t_{min}','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
legend([s1,s2], {'vs. prior','vs. posterior'});
ylabel(ax, 'Energy ratio');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
% ylim(ax,[0, 0.002]);
hold(ax,'off');

ax = f.ax(2);
hold(ax,'on'); 
scatter(ax, allrst_new(:,end), allrst_new(:,13), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
scatter(ax, allrst_new(:,end), allrst_new(:,14), msizehf, 'filled','bo');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
% plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'4-s','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Energy ratio');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
% ylim(ax,[0, 0.002]);
hold(ax,'off');

ax = f.ax(3);
hold(ax,'on'); 
scatter(ax, allrst_new(:,end), allrst_new(:,15), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
scatter(ax, allrst_new(:,end), allrst_new(:,16), msizehf, 'filled','bo');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
% plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'8-s','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Energy ratio');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
% ylim(ax,[0, 0.002]);
hold(ax,'off');

ax = f.ax(4);
hold(ax,'on'); 
scatter(ax, allrst_new(:,end), allrst_new(:,17), msizehf, 'filled','ro');  %, 'MarkerEdgeColor', 'w')
scatter(ax, allrst_new(:,end), allrst_new(:,18), msizehf, 'filled','bo');  %, 'MarkerEdgeColor', 'w')
% plot(ax, ax.XLim, [ttol1 ttol1], 'k--');
% plot(ax, ax.XLim, [ttol2 ttol2], 'b--');
text(ax,0.92,0.9,'16-s','FontSize',11,'unit','normalized','horizontalalignment','left',...
    'EdgeColor','k','Margin',2);
ylabel(ax, 'Energy ratio');
xlabel(ax, 'Time (day) since Mar. 1, 2003 (2003060)');
% ylim(ax,[0, 0.002]);
hold(ax,'off');

function [ax,hf] = plt_detections_inmap(ax,hftime,trange,crange,colnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ax] = plt_detections_inmap(ax,hftime,trange,colnum)
% This function is to plot the detections in the time range as input on the 
% map, color coded by time. The main part is similar to the first subplot of
% 'mig_linear_fit_LZB_addfam_autortm_v3.m'.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/10/25
% Last modified date:   2021/10/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

timecol = colnum(1:2);
%     indhf = find();
hf = hftime(hftime(:,timecol(1))==trange(1) & hftime(:,timecol(2))>=trange(2) & ...
            hftime(:,timecol(2))<=trange(3),:);


% subplot 1 of figure i
hold(ax,'on');
ax.FontSize = 9;
hfplt = sortrows(hf,-timecol(2));
xycol = colnum(3:end);
scatter(ax,hfplt(:,xycol(1)),hfplt(:,xycol(2)), 30, hfplt(:,timecol(2))/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
colormap(ax,'jet');
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.02];
%     c.TickLabels=[];
juldate = num2str(trange(1));
yr = str2double(juldate(1:4));
date = str2double(juldate(5:end));
a = jul2dat(yr,date);
mo = a(1);
if mo == 9
  mo = {' Sep. '};
elseif mo == 7
  mo = {' Jul. '};
else
  mo = {' Mar. '};
end
day = num2str(a(2));
yr = num2str(a(3));
c.Label.String = strcat({'Time (hr) on '}, day, mo, yr);
c.Label.FontSize = 11;
caxis(ax,[crange(1)/3600 crange(2)/3600]);

plot(ax,ax.XLim,[0 0],'k--');
plot(ax,[0 0],ax.YLim,'k--');
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');

% text(ax,0.85,0.1,'HF','FontSize',12,'unit','normalized');
% text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.04,0.1,'TWKB','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
medxhf = median(hfplt(:,1));
medyhf = median(hfplt(:,2));
%     xticks(f.ax(1),xran(1):5:xran(2));
%     yticks(f.ax(1),yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
%     text(f.ax(1),0.5,0.1,num2str(i),'FontSize',12,'unit','normalized');
text(ax,0.95,0.95,strcat(num2str(size(hfplt,1)),{' detections'}),'FontSize',10,...
  'unit','normalized','horizontalalignment','right');
text(ax,0.95,0.90,strcat({'in '},num2str(trange(3)-trange(2)),{' s'}),'FontSize',10,...
  'unit','normalized','horizontalalignment','right');
rate = sprintf('%.3f',size(hfplt,1)/(trange(3)-trange(2)));
text(ax,0.95,0.85,strcat({'rate: '},rate),'FontSize',10,...
  'unit','normalized','horizontalalignment','right');
hold(ax,'off');
    
% keyboard
    
    
    
    
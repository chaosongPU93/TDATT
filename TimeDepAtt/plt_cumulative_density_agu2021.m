function [f1] = plt_cumulative_density_agu2021(hfplt,xran,yran)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f1] = plt_cumulative_density(hfplt,lfplt,xran,yran)
% This function is to plot the cumulative density map of the hf and lf 
% catalog to see the activated region of the catalogs
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/28
% Last modified date:   2021/09/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[scrsz, res] = pixelperinch(1);

f1.fig=figure;
f1.fig.Renderer='Painters';
widin = 6;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f1.ax(isub) = subplot(nrow,ncol,isub);
end

set(f1.ax(1), 'position', [0.08 0.15 0.95-0.08 0.9-0.15]);

% xran = [-20 25];
% yran = [-20 20];

msizehf = 4;

% subplot 1 of figure i
% hfplt = hftime(hftime(:,13) > 2004*1000 & hftime(:,13) < 2005*1000, :);
% lfplt = lftime(lftime(:,13) > 2004*1000 & lftime(:,13) < 2005*1000, :);

ax = f1.ax(1);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
if ~isempty(hfplt)
  dxhf = 0.2;
  dyhf = 0.2;
  [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
  dumhf = density1d(density1d(:,3)>0, :);
  dumhf(dumhf(:,3)>1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
  dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
  dumhf(dumhf(:,3)==1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
  % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
  oldcmap = colormap(ax,'jet');
  % colormap(ax, flipud(oldcmap) );
  c=colorbar(ax,'SouthOutside');
  pos = ax.Position;
  c.Position = [pos(1), pos(2)-0.03, pos(3), 0.02];
  c.Label.String = strcat('log_{10}(N) of detections');
  c.Label.FontSize = 11;
% text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
  text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
  text(ax,0.1,0.93,strcat(num2str(length(hfplt(:,1))),{' detections'}),'FontSize',10,'unit','normalized',...
      'horizontalalignment','left');
end
ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):5:xran(2));
yticks(ax,yran(1):5:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);

%plot black box
bnd = [-10 -5;
       5 -5;
       5 5;
       -10 5;
       -10 -5;
      ];
plot(ax,bnd(:,1),bnd(:,2),'k-','linew',2);

% plot station names
dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');

stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4];
lzbtrio = stainfo(ind,:);
for i = 1: 3
    [dxlzb(i),dylzb(i)] = absloc2relaloc(str2num(lzbtrio(i,3)),str2num(lzbtrio(i,2)),...
                    -123.585000, 48.436667);
end

ind=[5,8,6];
pgctrio = stainfo(ind,:);
for i = 1: 3
    [dxpgc(i),dypgc(i)] = absloc2relaloc(str2num(pgctrio(i,3)),str2num(pgctrio(i,2)),...
                    -123.585000, 48.436667);
end

otherind = setdiff(1:12, union([3,11,4],[5,8,6]));
othersta = stainfo(otherind,:);
for i = 1: length(otherind)
    [dxother(i),dyother(i)] = absloc2relaloc(str2num(othersta(i,3)),str2num(othersta(i,2)),...
                    -123.585000, 48.436667);
end

%plot stations and names
for i=1:3
    % p1 is lzb trio stations
    p1 = plot(ax,dxlzb(i),dylzb(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','k');
    text(ax,dxlzb(i)+1,dylzb(i)+1,lzbtrio(i,1));
    % p2 is pgc trio stations
    p2 = plot(ax,dxpgc(i),dypgc(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','r');
    text(ax,dxpgc(i)+1,dypgc(i)+1,pgctrio(i,1)); 
end

for i = 1: length(otherind)
  % p3 is other unused stations
    p3 = plot(ax,dxother(i),dyother(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','w','clipping','on');
    text(ax,dxother(i)+1,dyother(i)+1,othersta(i,1),'clipping','on');
end

%plot fam names
loccont = [
           -123.585000 48.436667 36.8800;   % 002   
%            -123.549500 48.540833 38.5600;   % 243
%            -123.382333 48.574167 40.9800;   % 240
           ];
  nfampool = [
                '002';
%                 '243';
%                 '240';  % the most recent catalog of PGC is from fam 002, 243 and 240
              ];
         
[relacont(:,1), relacont(:,2)] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.585000, 48.436667);
scatter(ax,relacont(:,1), relacont(:,2),30,'k','filled');
for i = 1: size(loccont,1)
  text(ax,relacont(i,1)+1, relacont(i,2)+1, nfampool(i,:),'fontsize',8,'backgroundColor','w','Margin',0.5);
end

whitebox = [12 -8;
            49 -8;
            49 15;
            12 15;
            12 -8;
            ];
patch(ax,whitebox(:,1),whitebox(:,2),'w','edgecolor','none');          

hold(ax,'off');

set(f1.ax(2), 'position', [0.48 0.2 0.45 0.3]);
ax = f1.ax(2);
hold(ax,'on');
  dumhf = density1d(density1d(:,3)>0, :);
  dumhf(dumhf(:,3)>1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
  dumhf = sortrows(density1d(density1d(:,3)>0, :), 3);
  dumhf(dumhf(:,3)==1, :) = [];
  scatter(ax,dumhf(:,1),dumhf(:,2), msizehf, log10(dumhf(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
colormap(ax,'jet');
xlim(ax,[-10 5]);
ylim(ax,[-5 5]);
ax.Box = 'on';
ax.XTick = [];
ax.YTick = [];
axis(ax, 'equal');

keyboard













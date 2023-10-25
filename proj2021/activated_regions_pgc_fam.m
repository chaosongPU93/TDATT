% function activated_regions_pgc_fam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to use the density plot to lighten up the activated regions
% in HF and LF detected with PGC trio for any fam 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/06/14
% Last modified date:   2021/06/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all
clc

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forproj21');
isolfpath = strcat(workpath, '/Seisbasics/hypoinverse/forproj21');

% load files
fampool = ['002';
           '243';
           '240';
           '253';
           '036';
           '251';
           ];
nfam = size(fampool,1);

%%% for the locations of lfe families, by inversion from hypoinverse
%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/forproj21/chat00p.?
loccont = [
           -123.585000 48.436667 36.8800;   % 002   
           -123.549500 48.540833 38.5600;   % 243
           -123.382333 48.574167 40.9800;   % 240
           ];

for ifam = 3: nfam
    fam = fampool(ifam, :);
winlensechf = 4;

winlensec = 16;
loopoffmax = 4;
xcmaxAVEnmin = 0.45;

% SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
%                   '.',num2str(ccminlf));
% fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
% 
% fname = fname;
% hfmaplocall = load(fname);
% % 10 cols, format is:
% %   lon lat dep off12 off13 off12sec off13sec date strongest_arr_time center_of_win
% 
% SUFFIXlf = strcat('lf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
%                   '.',num2str(ccminlf));
% fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXlf);
% 
% fname = fname;
% lfmaplocall = load(fname);
% 
% % convert absolute loc to relative loc
% [ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
% lon0 = hfmaplocall(ind,1);
% lat0 = hfmaplocall(ind,2);
% if isempty(lon0)
%     lon0 = -123.5850;
%     lat0 = 48.4367;
% end
% 
% % relative coordinates to fam 043
% loc043 = [-123.772167 48.493000 35.5900;];
% [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),loc043(1),loc043(2));
% hftime = [dx dy hfmaplocall];     % now increases to 12 cols
% 
% [dx, dy] = absloc2relaloc(lfmaplocall(:,1),lfmaplocall(:,2),loc043(1),loc043(2));
% lftime = [dx dy lfmaplocall];     % now increases to 12 cols
% 
% % hfrelalocall = hftime;
% [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
% % hfrelalocall(:,1) = dx;
% % hfrelalocall(:,2) = dy;
% hftime = [dx dy hftime];    % now increases to 14 cols
% 
% % lfrelalocall = lftime;
% [dx, dy] = absloc2relaloc(lfmaplocall(:,1),lfmaplocall(:,2),lon0,lat0);
% % lfrelalocall(:,1) = dx;
% % lfrelalocall(:,2) = dy;
% lftime = [dx dy lftime];    % now increases to 14 cols    
% 
% % so now the first 2 cols are locations relative to fam 002, next 2 cols are relative to 043
% % 14 cols, format is:
% %   E(own) N(own) E(043) N(043) lon lat dep off12 off13 off12sec off13sec date 
% %   strongest_arr_time center_of_win
% 
% [dxloc0, dyloc0] = absloc2relaloc(lon0,lat0,loc043(1),loc043(2));
% 
% % sort according to day, sec
% hftime = sortrows(hftime, [14, 13]);
% lftime = sortrows(lftime, [14, 13]);

SUFFIXhf = strcat('up.hf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
                  '.',num2str(xcmaxAVEnmin));
SUFFIXlf = strcat('lf.time.',num2str(winlensechf),'_',num2str(winlensec),'.',num2str(loopoffmax),...
                  '.',num2str(xcmaxAVEnmin));
hffname = strcat(rstpath, '/evtloc.pgcfam.pj21.nodcut.',SUFFIXhf);
lffname = strcat(rstpath, '/evtloc.pgcfam.pj21.nodcut.',SUFFIXlf);

hffam = load(hffname);
lffam = load(lffname);

% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep + previous 50 cols
%%% UPDATED at 2021/06/23
%%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
%%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
%%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
%%%   for the rest, n=n+4;
%%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
%%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
%%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
%%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
%%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
%%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)


hftime = hffam(hffam(:, 13)==str2double(fam), :);
lftime = lffam(lffam(:, 13)==str2double(fam), :);

%%% load locations of isolated lf detections
fname = strcat(isolfpath, '/evtloc.offset_all_isolatedlf');
isolf = load(fname);
[dx, dy] = absloc2relaloc(isolf(:,1),isolf(:,2),-123.585000, 48.436667);
isolf = [dx dy isolf];
if isequal(fam, '002')
    isolf = isolf(1:2, :);
elseif isequal(fam, '243')
    isolf = isolf(3:10, :);
elseif isequal(fam, '240')
    isolf = isolf(11:12, :);
else
    isolf = [];
end

%convert to the same origin
[dxloc0, dyloc0] = absloc2relaloc(loccont(ifam,1),loccont(ifam,2),loccont(1,1),loccont(1,2));

%plotting range
if isequal(fam, '002')
    xran = [-35 15];
    yran = [-30 20];
elseif isequal(fam, '243')
    xran = [-30 20];
    yran = [-15 35];
elseif isequal(fam, '240')
    xran = [-20 30];
    yran = [-5 45];
else
    xran = [-35 15];
    yran = [-30 20];
end


% %% density in each ETS, 2003, 2004 and 2005
% f1.fig=figure;
% f1.fig.Renderer='Painters';
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 11;   % maximum height allowed is 11 inches
% set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% nrow = 3;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f1.ax(isub) = subplot(nrow,ncol,isub);
% end
% 
% % %%% reposition
% % set(f1.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
% % set(f1.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
% % set(f1.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
% % set(f1.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);
%     
% 
% msizehf = 4;
% msizelf = 10;
% 
% % subplot 1 of figure i
% hfplt = hftime(hftime(:,14) > 2003*1000 & hftime(:,14) < 2004*1000, :);
% lfplt = lftime(lftime(:,14) > 2003*1000 & lftime(:,14) < 2004*1000, :);
% 
% ax = f1.ax(1);
% hold(ax,'on');
% % plot(ax,[-100 100],[0 0],'k--');
% % plot(ax,[0 0],[-100 100],'k--');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxhf = 0.2;
% dyhf = 0.2;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% text(ax, 0.85, 0.93, '2003','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.5,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% % subplot 2 of figure i
% ax = f1.ax(2);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxlf = 0.5;
% dylf = 0.5;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% text(ax, 0.85, 0.93, '2003','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.5,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% 
% % subplot 3 of figure i
% hfplt = hftime(hftime(:,14) > 2004*1000 & hftime(:,14) < 2005*1000, :);
% lfplt = lftime(lftime(:,14) > 2004*1000 & lftime(:,14) < 2005*1000, :);
% 
% ax = f1.ax(3);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxhf = 0.2;
% dyhf = 0.2;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.5,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% % subplot 4 of figure i
% ax = f1.ax(4);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxlf = 0.5;
% dylf = 0.5;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% text(ax, 0.85, 0.93, '2004','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.5,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% 
% % subplot 5 of figure i
% hfplt = hftime(hftime(:,14) > 2005*1000, :);
% lfplt = lftime(lftime(:,14) > 2005*1000, :);
% 
% ax = f1.ax(5);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxhf = 0.2;
% dyhf = 0.2;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% text(ax, 0.85, 0.93, '2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'c','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.5,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% % subplot 6 of figure i
% ax = f1.ax(6);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% 
% % create a density matrix to store the number of detections in each small grid
% dxlf = 0.5;
% dylf = 0.5;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% text(ax, 0.85, 0.93, '2005','FontSize',12,'unit','normalized','horizontalalignment','center',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'d','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.5,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% print(f1.fig,'-dpdf',strcat(rstpath,'/density_eachETS_',fam,'.pdf'));
% 
% 
% %% 2003 + 2004 + 2005 combined
% f2.fig=figure;
% f2.fig.Renderer='Painters';
% widin = 8;  % maximum width allowed is 8.5 inches
% htin = 9;   % maximum height allowed is 11 inches
% set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
% nrow = 2;
% ncol = 2;
% for isub = 1:nrow*ncol
%     f2.ax(isub) = subplot(nrow,ncol,isub);
% end
% 
% % %%% reposition
% % set(f2.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
% % set(f2.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
% % set(f2.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
% % set(f2.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);
% 
% % xran = [-35 15];
% % yran = [-30 20];
% 
% msizehf = 4;
% msizelf = 10;
% 
% % subplot 1 of figure i
% hfplt = hftime(hftime(:,14) > 2003*1000 & hftime(:,14) < 2004*1000, :);
% lfplt = lftime(lftime(:,14) > 2003*1000 & lftime(:,14) < 2004*1000, :);
% 
% ax = f2.ax(1);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% %     ind = find(density1d>=0);
% %     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% %     aaa = density2d';
% %     bbb = aaa(aaa>=1);
% %     ccc = log(bbb);
% %     imagesc(ax,xran+dx, yran+dy, ccc);
% % create a density matrix to store the number of detections in each small grid
% dxhf = 0.2;
% dyhf = 0.2;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
% text(ax, 0.95, 0.93, '2003','FontSize',12,'unit','normalized','horizontalalignment','right',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.4,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% % subplot 2 of figure i
% ax = f2.ax(2);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxlf = 0.5;
% dylf = 0.5;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
% % scatter(ax,isolf(:,1),isolf(:,2),15,'k');
% text(ax, 0.95, 0.93, '2003','FontSize',12,'unit','normalized','horizontalalignment','right',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.4,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):10:xran(2));
% yticks(ax,yran(1):10:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% 
% % subplot 3 of figure i
% hfplt = hftime(hftime(:,14) > 2004*1000, :);
% lfplt = lftime(lftime(:,14) > 2004*1000, :);
% 
% ax = f2.ax(3);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% %     ind = find(density1d>=0);
% %     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% %     aaa = density2d';
% %     bbb = aaa(aaa>=1);
% %     ccc = log(bbb);
% %     imagesc(ax,xran+dx, yran+dy, ccc);
% % create a density matrix to store the number of detections in each small grid
% dxhf = 0.2;
% dyhf = 0.2;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
% text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.4,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% % subplot 4 of figure i
% ax = f2.ax(4);
% hold(ax,'on');
% scatter(ax,0,0,15,'k','filled')
% plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
% plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
% ax.FontSize = 9;
% % create a density matrix to store the number of detections in each small grid
% dxlf = 0.5;
% dylf = 0.5;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
% dum = sortrows(density1d(density1d(:,3)>0, :), 3);
% dum(dum(:,3)==1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% % imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
% oldcmap = colormap(ax,'jet');
% % colormap(ax, flipud(oldcmap) );
% c=colorbar(ax,'SouthOutside');
% pos = ax.Position;
% c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
% c.Label.String = strcat('log_{10}(N) of detections');
% c.Label.FontSize = 8;
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
%     xran,yran,dxlf,dylf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
% % scatter(ax,isolf(:,1),isolf(:,2),15,'k','linew',1.5);
% text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
%      'EdgeColor','k','Margin',2);
% text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
% text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
% text(ax,0.4,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
%      'horizontalalignment','center');
% ax.Box = 'on';
% grid(ax, 'on');
% axis(ax, 'equal');
% ax.GridLineStyle = '--';
% ax.XAxisLocation = 'top';
% xlim(ax,xran);
% ylim(ax,yran);
% xticks(ax,xran(1):5:xran(2));
% yticks(ax,yran(1):5:yran(2));
% xlabel(ax,'E (km)','fontsize',11);
% ylabel(ax,'N (km)','fontsize',11);
% hold(ax,'off');
% 
% print(f2.fig,'-dpdf',strcat(rstpath,'/density_combinedETS_',fam,'.pdf'));



%% 2004 + 2005 combined, excluding 2003
%%% as all excellent LF detections where the time lag from min to max is very short compared to full
%%% duration
f2.fig=figure;
f2.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 9;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 2;
ncol = 1;
for isub = 1:nrow*ncol
    f2.ax(isub) = subplot(nrow,ncol,isub);
end

% %%% reposition
% set(f2.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
% set(f2.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
% set(f2.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
% set(f2.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);

% xran = [-35 15];
% yran = [-30 20];

msizehf = 4;
msizelf = 10;


% subplot 1 of figure i
hfplt = hftime(hftime(:,14) > 2004*1000, :);
lfplt = lftime(lftime(:,14) > 2004*1000, :);

ax = f2.ax(1);
hold(ax,'on');
scatter(ax,0,0,15,'k','filled')
plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
ax.FontSize = 9;
%     ind = find(density1d>=0);
%     scatter(ax,xloc(ind),yloc(ind), 4, log(density1d(ind)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
%     aaa = density2d';
%     bbb = aaa(aaa>=1);
%     ccc = log(bbb);
%     imagesc(ax,xran+dx, yran+dy, ccc);
% create a density matrix to store the number of detections in each small grid
dxhf = 0.2;
dyhf = 0.2;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,1),hfplt(:,2),...
    xran,yran,dxhf,dyhf);
dum = density1d(density1d(:,3)>0, :);
dum(dum(:,3)>1, :) = [];
scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(density1d(density1d(:,3)>0, :), 3);
dum(dum(:,3)==1, :) = [];
scatter(ax,dum(:,1),dum(:,2), msizehf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 8;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
    xran,yran,dxhf,dyhf);
dum = density1d(density1d(:,3)>0, :);
dum(dum(:,3)>1, :) = [];
scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
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
hold(ax,'off');

% subplot 2 of figure i
ax = f2.ax(2);
hold(ax,'on');
scatter(ax,0,0,15,'k','filled')
plot(ax,[-100 100],[dyloc0 dyloc0],'k--');
plot(ax,[dxloc0 dxloc0],[-100 100],'k--');
ax.FontSize = 9;
% create a density matrix to store the number of detections in each small grid
dxlf = 0.5;
dylf = 0.5;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(lfplt(:,1),lfplt(:,2),...
    xran,yran,dxlf,dylf);
dum = density1d(density1d(:,3)>0, :);
dum(dum(:,3)>1, :) = [];
scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)),'s','linew',0.2);  %, 'MarkerEdgeColor', 'w')
dum = sortrows(density1d(density1d(:,3)>0, :), 3);
dum(dum(:,3)==1, :) = [];
scatter(ax,dum(:,1),dum(:,2), msizelf, log10(dum(:,3)), 'filled','s');  %, 'MarkerEdgeColor', 'w')
% imagesc(ax,[xran(1)+0.5*dx xran(end)-0.5*dx], [yran(1)+0.5*dy yran(end)-0.5*dy], density2d');
oldcmap = colormap(ax,'jet');
% colormap(ax, flipud(oldcmap) );
c=colorbar(ax,'SouthOutside');
pos = ax.Position;
c.Position = [pos(1), pos(2)-0.018, pos(3), 0.015];
c.Label.String = strcat('log_{10}(N) of detections');
c.Label.FontSize = 8;
[density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
    xran,yran,dxlf,dylf);
dum = density1d(density1d(:,3)>0, :);
dum(dum(:,3)>1, :) = [];
scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
% scatter(ax,isolf(:,1),isolf(:,2),15,'k','linew',1.5);
text(ax, 0.95, 0.93, '2004+2005','FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'LF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'b','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(lfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');
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
hold(ax,'off');

print(f2.fig,'-dpdf',strcat(rstpath,'/density_combinedETS0405_',fam,'.pdf'));


end
















function f = pltevtloc(f,regflag,evtall,tmrall,slab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function plot the spatial distribution of the regular events 
% and tremor, color-coded by different magnitudes, along with the depth 
% variation by a pseudo 3-D view through latitude and longitude 
% 
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/18
% Last modified date:   2020/02/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtain the depth grid file for contour
lat = 31:0.01:36;
lon = 131:0.01:138;
[longrd, latgrd] = meshgrid(lon,lat);
depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');

F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
deptmr = F(tmrall(:,6),tmrall(:,7));

% sort the event according to the magnitude
evtall = sortrows(evtall,11);

% subplot 1, horizontal distribution
hold(f.ax(1),'on');
f.ax(1).Box = 'on';
grid(f.ax(1), 'on');
if regflag == 1  % means western shikoku
    plot(f.ax(1),[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
         [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f.ax(1),[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
% axis(f.ax(1), 'equal');
if regflag == 1
    xlim(f.ax(1),[131 135]);
    ylim(f.ax(1),[32.3 35.2]);
elseif regflag == 2
    xlim(f.ax(1),[134.5 137.2]);
    ylim(f.ax(1),[33 35.5]);
end
contour(f.ax(1),longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',250,'linew',1,...
    'linec',[0.8 0.8 0.8]);
coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(f.ax(1),coast(:,1),coast(:,2),'black','linew',0.5);
scatter(f.ax(1),tmrall(:,6),tmrall(:,7),2,[0.4 0.4 0.4],'filled','o');
scatter(f.ax(1),evtall(:,8),evtall(:,9),8,evtall(:,11),'filled','o');
colormap(f.ax(1),'jet');
c=colorbar(f.ax(1),'northoutside');
pos = f.ax(1).Position;
c.Label.String = 'Magnitude';
c.Label.FontSize = 10;
caxis(f.ax(1),[min(evtall(:,11)) max(evtall(:,11))]);
xlabel(f.ax(1),{'Longitude ({\circ})'});
ylabel(f.ax(1),{'Latitude ({\circ})'});
hold(f.ax(1),'off');

% subplot 2, vertical distribution along latitude
pos1 = f.ax(1).Position;
pos2 = f.ax(2).Position;
f.ax(2).Position = [pos2(1) pos1(2) pos2(3) pos1(4)];
hold(f.ax(2),'on');
f.ax(2).Box = 'on';
grid(f.ax(2), 'on');
if regflag == 1
    xlim(f.ax(2),[10 100]);
    ylim(f.ax(2),[32.3 35.2]);
    xticks(f.ax(2),10:10:100);
elseif regflag == 2
    xlim(f.ax(2),[10 80]);
    ylim(f.ax(2),[33 35.5]);
    xticks(f.ax(2),10:10:80);
end
scatter(f.ax(2),deptmr,tmrall(:,7),2,[0.4 0.4 0.4],'filled','o');
scatter(f.ax(2),evtall(:,10),evtall(:,9),8,evtall(:,11),'filled','o');
colormap(f.ax(2),'jet');
% c=colorbar(f.ax(2),'northoutside');
caxis(f.ax(2),[min(evtall(:,11)) max(evtall(:,11))]);
xlabel(f.ax(2),{'Depth (km)'});
hold(f.ax(2),'off');

% subplot 3, vertical distribution along longitude
pos2 = f.ax(2).Position;
pos3 = f.ax(3).Position;
f.ax(3).Position = [pos1(1) pos3(2) pos1(3) pos2(3)];
hold(f.ax(3),'on');
f.ax(3).Box = 'on';
grid(f.ax(3), 'on');
if regflag == 1
    xlim(f.ax(3),[131 135]);
    ylim(f.ax(3),[10 100]);
    yticks(f.ax(3),10:10:100);
elseif regflag == 2
    xlim(f.ax(3),[134.5 137.2]);
    ylim(f.ax(3),[10 80]);
    yticks(f.ax(3),10:10:80);
end
set(f.ax(3), 'YDir','reverse');
scatter(f.ax(3),tmrall(:,6),deptmr,2,[0.4 0.4 0.4],'filled','o');
scatter(f.ax(3),evtall(:,8),evtall(:,10),8,evtall(:,11),'filled','o');
colormap(f.ax(3),'jet');
% c=colorbar(f.ax(3),'northoutside');
caxis(f.ax(3),[min(evtall(:,11)) max(evtall(:,11))]);
ylabel(f.ax(3),{'Depth (km)'});
hold(f.ax(3),'off');

% subplot 4, horizontal distribution
hold(f.ax(4),'on');
f.ax(4).Box = 'on';
grid(f.ax(4), 'on');
if regflag == 1  % means western shikoku
    plot(f.ax(4),[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
         [0.6 0.6 0.6],'linew',2);
elseif regflag == 2  % means kii pen
    plot(f.ax(4),[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
        [0.6 0.6 0.6],'linew',2);
end
% axis(f.ax(4), 'equal');
if regflag == 1
    xlim(f.ax(4),[131 135]);
    ylim(f.ax(4),[32.3 35.2]);
elseif regflag == 2
    xlim(f.ax(4),[134.5 137.2]);
    ylim(f.ax(4),[33 35.5]);
end
contour(f.ax(4),longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',250,'linew',1,...
    'linec',[0.8 0.8 0.8]);
coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(f.ax(4),coast(:,1),coast(:,2),'black','linew',0.5);
scatter(f.ax(4),evtall(:,8),evtall(:,9),8,evtall(:,11),'filled','o');
scatter(f.ax(4),tmrall(:,6),tmrall(:,7),2,[0.4 0.4 0.4],'filled','o');
colormap(f.ax(4),'jet');
c=colorbar(f.ax(4),'northoutside');
pos = f.ax(4).Position;
c.Label.String = 'Magnitude';
c.Label.FontSize = 10;
caxis(f.ax(4),[min(evtall(:,11)) max(evtall(:,11))]);
xlabel(f.ax(4),{'Longitude ({\circ})'});
ylabel(f.ax(4),{'Latitude ({\circ})'});
hold(f.ax(4),'off');

% keyboard











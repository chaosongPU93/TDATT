function ax = pltevttmrloc2(ax,regflag,slab,evtobj,tmrobj,minmag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to plot the location of tremor detections and regular
% events during the same time period, color-code the events by magnitude
% 
%   INPUT:
%       ax:     axes to plot on
%       regflag: flag to indicate the region, 1 for western shikoku, 2 for kii
%                penninsula
%       slab:   slab grid file
%       evtobj:  object regular event 
%       tmrobj:  object tremor 
%       magtol:     magnitude threshold
%   OUTPUT:
%       ax:     return the same axes as input
%       
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/26
% Last modified date:   2020/02/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% obtain the depth grid file for contour
lat = 31:0.01:36;
lon = 131:0.01:138;
[longrd, latgrd] = meshgrid(lon,lat);
depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');

% apply a minimum magnitude threshold if needed
evtobj = evtobj(evtobj(:,11)>=minmag, :);

% sort the event according to the magnitude
evtobj = sortrows(evtobj,11);

%
hold(ax,'on');
if regflag == 1  % means western shikoku
    plot(ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',[0.6 0.6 0.6],...
         'linew',2);
elseif regflag == 2  % means kii pen
    plot(ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',[0.6 0.6 0.6],...
         'linew',2);
end
contour(ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
        'linec',[0.8 0.8 0.8]);
coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(ax,coast(:,1),coast(:,2),'black','LineWidth',0.5);
scatter(ax,tmrobj(:,6),tmrobj(:,7),2,[0.4 0.4 0.4],'filled','o');
scatter(ax,evtobj(:,8),evtobj(:,9),8,evtobj(:,11),'filled','o');
colormap(ax,'jet');
c=colorbar(ax,'northoutside');
c.Label.String = 'Magnitude';
c.Label.FontSize = 10;
% caxis(ax,[min(evtobj(:,11)) max(evtobj(:,11))]);
if minmag > -5
    caxis(ax,[minmag 5]);
else
    caxis(ax,[-0.5 5]);
end
ax.Box = 'on';
grid(ax, 'on');
% axis(ax, 'equal');

if regflag == 1
    xlim(ax,[131 135]);
    ylim(ax,[32.3 35.2]);
elseif regflag == 2
    xlim(ax,[134.5 137.2]);
    ylim(ax,[33 35.5]);
end
hold(ax,'off');
% % The data aspect ratio is the relative length of the data units along the x-axis, y-axis, and z-axis.
% lat0 = 0.5*(ax.YLim(1)+ax.YLim(2));
% set(ax,'DataAspectRatio',[1/cosd(lat0) 1 1])


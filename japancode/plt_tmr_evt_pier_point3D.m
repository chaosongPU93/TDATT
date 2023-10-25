function ax = plt_tmr_evt_pier_point3D(ax,regflag,slab,tmr,pierpt,stasel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to plot the 3D spatial distribution of tremors and the
% piercing points of the events at the slab interface. The tremor location
% is not variant regardless of the stations. In contrast, the piercing points
% depend on the station. So this function accepts the station query.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/09
% Last modified date:   2020/04/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% % get the scrsz in pixels and number of pixels per inch of monitor 1
% [scrsz, res] = pixelperinch(1);
% 
% workpath = '/home/data2/chaosong/shikoku_kii';
% figpath = strcat(workpath,'/figs');

if regflag == 1
    prefix = 'shikoku';
elseif regflag == 2
    prefix = 'kii';
end


%% interpolate the slab interface to slightly denser grid
if regflag == 1
    lat = 32.5:0.02:35.5;
    lon = 131.5:0.02:135;
elseif regflag == 2
    lat = 32.5:0.02:35;
    lon = 134:0.02:137.5;
end
[longrd, latgrd] = meshgrid(lon,lat);

% create a interpolant object function for multiple usages
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
depgrd = F(longrd,latgrd);
% depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'cubic');

% % for plotting purpose
% slabspa = [reshape(longrd,[],1) reshape(latgrd,[],1) reshape(depgrdspa,[],1)];
% [slabspac(:,1),slabspac(:,2)] = absloc2relaloc(slabspa(:,1),slabspa(:,2),lon0,lat0);
% slabspac(:,3) = slabspa(:,3);
% longrdspac = reshape(slabspac(:,1),size(longrd));
% latgrdspac = reshape(slabspac(:,2),size(latgrd));
% depgrdspac = reshape(slabspac(:,3),size(depgrdspa));


%% interpolate the tremors on the slab interface
tmr(:,8) = F(tmr(:,6),tmr(:,7));


%%
% nsta = size(pierpt,1);
% nevt = size(pierpt,2);
% 
% i = stanum;
% pierptevt(:,:) = pierpt(i,:,:);

% ax=gca;
hold(ax,'on');
ax.Box = 'on';
grid(ax,'on');
surf(ax,longrd,latgrd,depgrd,'edgecolor','none');
p1=scatter3(ax,tmr(:,6),tmr(:,7),tmr(:,8),5,[0.6 0.6 0.6],'filled','o');
p2=scatter3(ax,pierpt(:,1),pierpt(:,2),pierpt(:,3),5,[1 0.6 0],'filled','o');
scatter3(ax,stasel.lo,stasel.la,0,60,'y','filled','^','markeredgec','k');
legend(ax,[p1,p2],{'Tremors','Piercing points of events'},'fontsize',10);
set(ax,'ZDir','reverse');
colorbar;
view(ax,185,15);
xlabel(ax,'Longitude ^o');
ylabel(ax,'Latitude ^o');
zlabel(ax,'Depth (km)');

% title(ax,strcat(prefix,'---',stasel.nm));

hold(ax,'off');
  
% keyboard








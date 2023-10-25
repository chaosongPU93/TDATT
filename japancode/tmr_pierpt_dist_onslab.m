function [ind,eucdist,ind2,eucdist2,numtmr,kdtreetmr,f1,f2,f3,f4,f5]=...
    tmr_pierpt_dist_onslab(regflag,plotflag,slab,tmr,evt,pierpt,sta,range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to find the distance from the tremor sources to the
% piercing points of events calculated on the station base. It applies the
% kdtree model to all tremors, and then use knnsearch to get the k-nearest
% tremor sources to each events, ie, it can find if how close is each
% individual events to the tremor cluster, if this closest distance is larger
% than a threshold, then it means this event could never share the same path
% as any tremor, assuming that the tremor cluster indicates a slipping zone.
% Additionally, it uses a rangesearch, which searches all tremor sources that
% are within this range to each individual event, assuming tremor sources
% outside this range would not have slipping effect on this event.
%
%
%
% kdtree:
% 1.Kd-trees divide your data into nodes with at most BucketSize (default
%   is 50) points per node, based on coordinates (as opposed to categories).
%
%
% knn search:
% 1.Given a set X of n points and a distance function, k-nearest neighbor (kNN)
%   search lets you find the k closest points in X to a query point or set of
%   points Y.
% 2.In contrast, for a positive real value r, rangesearch finds all points in X
%   that are within a distance r of each point in Y. This fixed-radius search is
%   closely related to kNN search, as it supports the same distance metrics and
%   search classes, and uses the same search algorithms.
% 3.knnsearch does the following: Determines the node to which the query point
%   belongs. Finds the closest k points within that node and its distance to the
%   query point. Chooses all other nodes having any area that is within the same
%   distance, in any direction, from the query point to the kth closest point.
%   Searches nodes within that range for any points closer to the query point.
%
%
% range search:
% 1. Find all neighbors within specified distance using searcher object.
% 2. Basically the same syntax as knnsearch, Idx = rangesearch(Mdl,Y,r) searches
%   for all neighbors (i.e., points, rows, or observations) in Mdl.X within radius
%   r of each point (i.e., row or observation) in the query data Y using an
%   exhaustive search or a Kd-tree. rangesearch returns Idx, which is a column
%   vector of the indices of Mdl.X within r units.
%
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/16
% Last modified date:   2020/04/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = '/home/data2/chaosong/shikoku_kii';
datapath = strcat(workpath,'/matsave');

%% create a kdtree model of slab interface in cartesian coordinate
if regflag == 1
    prefix = 'shikoku';
    lat = 32.5:0.01:35.5;
    lon = 131.5:0.01:135;
elseif regflag == 2
    prefix = 'kii';
    lat = 32.5:0.01:35;
    lon = 134:0.01:137.5;
end
lat0 = 0.5*(lat(1)+lat(end));
lon0 = 0.5*(lon(1)+lon(end));
[tmrr(:,1),tmrr(:,2)] = absloc2relaloc(tmr(:,6),tmr(:,7),lon0,lat0);

% create a interpolant object function for multiple usages
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
tmrr(:,3) = F(tmr(:,6),tmr(:,7));

% create a kd tree for all tremors
kdtreetmr = KDTreeSearcher(tmrr,'Distance','euclidean','BucketSize',50);

%%% also convert to the piercing points from geographical to cartesian coord
[pierptr(:,1),pierptr(:,2)] = absloc2relaloc(pierpt(:,1),pierpt(:,2),lon0,lat0);
pierptr(:,3) = pierpt(:,3);


%% knnsearch, find the nearest tremor source to each piercing point

% ind, the length of ind should be the same as event number
% eucdist, the euclidean distance between each event to the closest tremor
[ind,eucdist] = knnsearch(kdtreetmr,pierptr,'K',1);

% the object closest tremor
tmrobj = tmr(ind, :);  % this index applies to tremor cluster


%% rangesearch, find the all tremor sources within the dist range to each piercing point
% i.e. within this radius, we think the ray path for all tremors is the same as the event at center

% range = 5;      % max dist, [km]
% the result is a cell array
[ind2,eucdist2] = rangesearch(kdtreetmr,pierptr,range);
numtmr = zeros(size(eucdist2));
mindist = nan(size(eucdist2));
for i = 1: size(eucdist2,1)
    numtmr(i) = size(eucdist2{i},2);
    if numtmr(i) ~= 0
        mindist(i) = min(eucdist2{i});
    end
end

if plotflag
    %% plot the distance distribution of events to their closest tremor source
    f1.fig=figure;
    f1.fig.Renderer='Painters';
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    set(f1.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    f1.ax = gca;
    hold(f1.ax,'on');
    f1.ax.Box = 'on';
    grid(f1.ax, 'on');
    if regflag == 1  % means western shikoku
        plot(f1.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    elseif regflag == 2  % means kii pen
        plot(f1.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    end
    axis(f1.ax, 'equal');
    if regflag == 1
        xlim(f1.ax,[131.5 133.5]);
        ylim(f1.ax,[32.5 34]);
        %     xticks(f.ax,131.8:0.25:134.6);
        %     yticks(f.ax,32.8:0.25:34.9);
        %     f.ax.XMinorGrid = 'on';
        %     f.ax.XMinorTick = 'on';
        %     f.ax.YMinorGrid = 'on';
        %     f.ax.YMinorTick = 'on';
    elseif regflag == 2
        xlim(f1.ax,[134.5 136.4]);
        ylim(f1.ax,[33.1 34.8]);
        %     xticks(f.ax,134.5:0.25:137);
        %     yticks(f.ax,33.1:0.25:35.4);
        %     f.ax.XMinorGrid = 'on';
        %     f.ax.XMinorTick = 'on';
        %     f.ax.YMinorGrid = 'on';
        %     f.ax.YMinorTick = 'on';
    end
    
    % lat = 31:0.01:36;
    % lon = 131:0.01:138;
    % [longrd, latgrd] = meshgrid(lon,lat);
    % depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');
    % contour(f.ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',300,'linew',1,...
    %         'linec',[0.8 0.8 0.8]);
    
    % scatter(f.ax,evt(eucdist>5,8),evt(eucdist>5,9),2,[0.5 0.5 0.5],'filled','o');
    %
    % scatter(f.ax,evt(eucdist<=5,8),evt(eucdist<=5,9),2,eucdist(eucdist<=5),'filled','o');
    
    scatter(f1.ax,pierpt(:,1),pierpt(:,2),2,eucdist,'filled','o');% colored by dist from closest tremor
%     scatter(f1.ax,evt(:,8),evt(:,9),2,eucdist,'filled','o');
    
    colormap(f1.ax,'jet');
    c=colorbar(f1.ax);
    c.Label.String = strcat('Min. dist to tremors via knnsearch (km)');
    c.Label.FontSize = 12;
    caxis(f1.ax,[0 5]);
    
    scatter(f1.ax,tmr(:,6),tmr(:,7),2,[0.8 0.8 0.8],'filled','o');
    
    scatter(f1.ax,tmrobj(:,6),tmrobj(:,7),2,[0.2 0.2 0.2],'filled','o');
    
    scatter(f1.ax,sta.lo, sta.la,80,'y','filled','^','markeredgec','k');
    
    text(f1.ax,sta.lo-0.2, sta.la+0.1, sta.nm, 'fontsize',12, 'color','k');
    
    % lat = 33.1:0.01:34.8;
    % lon = 134.5:0.01:136.4;
    % [longrd, latgrd] = meshgrid(lon,lat);
    % eucdistgrd = griddata(evt(:,8),evt(:,9),eucdist,longrd,latgrd, 'cubic');
    
    % contour(f.ax,longrd,latgrd,eucdistgrd,1:2:5,'showtext','on','LabelSpacing',400,'linew',0.5,...
    %         'linec',[0.5 0.5 0.5]);
    
    
    %% plot the histogram of the distance distribution
    f2.fig=figure;
    f2.fig.Renderer='Painters';
    widin = 4;  % maximum width allowed is 8.5 inches
    htin = 5;   % maximum height allowed is 11 inches
    set(f2.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    f2.ax = gca;
    hold(f2.ax, 'on');
    f2.ax.Box = 'on';
    grid(f2.ax, 'on');
    f2.ax.GridLineStyle = '--';
    patarea = [0 0;
        2.5 0;
        2.5 1;
        0 1;
        0 0];
    patch(f2.ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
    h = histogram(f2.ax,eucdist','binwidth',0.5,'normalization','probability','facecolor','k',...
        'facea',0.5);
    for i = 1: 10
        text(f2.ax,h.BinEdges(i),h.Values(i)+0.005,sprintf('%.2f',h.Values(i)));
    end
    perc = sum(eucdist<=2.5)/size(eucdist,1);
    text(f2.ax,0.25,0.9,sprintf('%.2f',perc),'fontsize',12,'unit','normalized');
    text(f2.ax,0.45,0.95,num2str(size(eucdist,1)),'fontsize',13,'unit','normalized');
    xlim(f2.ax,[0 5]);
    ylim(f2.ax,[0 0.26]);
    xlabel(f2.ax,'Distance (km)');
    ylabel(f2.ax,'Percentage');
    hold(f2.ax, 'off');
    
    
    %% plot the number of tremors that within search range for each event
    f3.fig=figure;
    f3.fig.Renderer='Painters';
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    set(f3.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    f3.ax = gca;
    hold(f3.ax,'on');
    f3.ax.Box = 'on';
    grid(f3.ax, 'on');
    if regflag == 1  % means western shikoku
        plot(f3.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    elseif regflag == 2  % means kii pen
        plot(f3.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    end
    axis(f3.ax, 'equal');
    if regflag == 1
        xlim(f3.ax,[131.5 133.5]);
        ylim(f3.ax,[32.5 34]);
        
    elseif regflag == 2
        xlim(f3.ax,[134.5 136.4]);
        ylim(f3.ax,[33.1 34.8]);
        
    end
    scatter(f3.ax,tmr(:,6),tmr(:,7),2,[0.8 0.8 0.8],'filled','o');
    
    
    scatter(f3.ax,pierpt(numtmr==0,1),pierpt(numtmr==0,2),2,[0.5 0.5 0.5],'filled','o');
    
    scatter(f3.ax,pierpt(numtmr>0,1),pierpt(numtmr>0,2),2,log(numtmr(numtmr>0)),'filled','o');
    
%     scatter(f3.ax,evt(numtmr==0,8),evt(numtmr==0,9),2,[0.5 0.5 0.5],'filled','o');
%     
%     scatter(f3.ax,evt(numtmr>0,8),evt(numtmr>0,9),2,log(numtmr(numtmr>0)),'filled','o');
    
    colormap(f3.ax,'jet');
    c=colorbar(f3.ax);
    c.Label.String = strcat('log_{10}(N) of tremors within 2.5 km of each piercing point');
    c.Label.FontSize = 12;
    % caxis(f3.ax,[1 5]);
    
    scatter(f3.ax,sta.lo, sta.la,80,'y','filled','^','markeredgec','k');
    
    text(f3.ax,sta.lo-0.2, sta.la+0.1, sta.nm, 'fontsize',12, 'color','k');
    
    
    %% plot the min. distance between tremor and each event from range search
    f4.fig=figure;
    f4.fig.Renderer='Painters';
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 8;   % maximum height allowed is 11 inches
    set(f4.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    f4.ax = gca;
    hold(f4.ax,'on');
    f4.ax.Box = 'on';
    grid(f4.ax, 'on');
    if regflag == 1  % means western shikoku
        plot(f4.ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    elseif regflag == 2  % means kii pen
        plot(f4.ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',...
            [0.6 0.6 0.6],'linew',2);
    end
    axis(f4.ax, 'equal');
    if regflag == 1
        xlim(f4.ax,[131.5 133.5]);
        ylim(f4.ax,[32.5 34]);
        
    elseif regflag == 2
        xlim(f4.ax,[134.5 136.4]);
        ylim(f4.ax,[33.1 34.8]);
        
    end
    scatter(f4.ax,tmr(:,6),tmr(:,7),2,[0.8 0.8 0.8],'filled','o');
    
    scatter(f4.ax,pierpt(numtmr==0,1),pierpt(numtmr==0,2),2,[0.5 0.5 0.5],'filled','o');
    
    scatter(f4.ax,pierpt(numtmr>0,1),pierpt(numtmr>0,2),2,mindist(numtmr>0),'filled','o');

%     scatter(f4.ax,evt(numtmr==0,8),evt(numtmr==0,9),2,[0.5 0.5 0.5],'filled','o');
%     
%     scatter(f4.ax,evt(numtmr>0,8),evt(numtmr>0,9),2,mindist(numtmr>0),'filled','o');
    
    colormap(f4.ax,'jet');
    c=colorbar(f4.ax);
    c.Label.String = strcat('Min. dist to tremors via within 2.5 km range (km)');
    c.Label.FontSize = 12;
    caxis(f4.ax,[0 5]);
    
    scatter(f4.ax,sta.lo, sta.la,80,'y','filled','^','markeredgec','k');
    
    text(f4.ax,sta.lo-0.2, sta.la+0.1, sta.nm, 'fontsize',12, 'color','k');
    
    
    %% plot the histogram of the number of tremors in range 
    f5.fig=figure;
    f5.fig.Renderer='Painters';
    widin = 4;  % maximum width allowed is 8.5 inches
    htin = 5;   % maximum height allowed is 11 inches
    set(f5.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    f5.ax = gca;
    hold(f5.ax, 'on');
    f5.ax.Box = 'on';
    grid(f5.ax, 'on');
    f5.ax.GridLineStyle = '--';
    % patarea = [0 0;
    %            2.5 0;
    %            2.5 1;
    %            0 1;
    %            0 0];
    % patch(f5.ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
    h2 = histogram(f5.ax,log(numtmr(numtmr>0)),'binwidth',1,'normalization','count',...
        'facecolor','k','facea',0.5);
    perc = sum(numtmr>0)/size(eucdist2,1);
    for i = 1: length(h2.Values)
        text(f5.ax,h2.BinEdges(i)+0.1,h2.Values(i),sprintf('%.2f',h2.Values(i)/size(eucdist2,1)));
    end
    text(f5.ax,0.2,0.95,num2str(size(eucdist2,1)),'fontsize',13,'unit','normalized');
    text(f5.ax,0.4,0.95,sprintf('%.2f',perc),'fontsize',13,'unit','normalized');
    % xlim(f5.ax,[0 5]);
    % ylim(f5.ax,[0 0.22]);
    xlabel(f5.ax,'log_{10}(N) of tremors');
    ylabel(f5.ax,'Count');
    hold(f5.ax, 'off');
    
else
    f1=[];
    f2=[];
    f3=[];
    f4=[];
    f5=[];
end

% keyboard


















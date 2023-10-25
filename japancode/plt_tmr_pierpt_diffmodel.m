function [f1,f2,diffdist]=plt_tmr_pierpt_diffmodel(regflag,plotflag,pierptref,pierptcomp,tmr,evt,sta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to plot the difference in piercing points
% of the same events in km between 2 velocity models, this makes more sense
% comparing to 'plt_tmr_pierpt_dist_diff.m' since the piercing points of same
% event under different models not only differ themselves, but also link to
% different closest tremor!
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/05/07
% Last modified date:   2020/05/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

%% calculate the distance between piercing points in km
if regflag == 1
    prefix = 'shikoku';
    lat = 32.5:0.1:35.5;
    lon = 131.5:0.1:135;
elseif regflag == 2
    prefix = 'kii';
    lat = 32.5:0.1:35;
    lon = 134:0.1:137.5;
end

lat0 = 0.5*(lat(1)+lat(end));
lon0 = 0.5*(lon(1)+lon(end));

[pierptref(:,4),pierptref(:,5)] = absloc2relaloc(pierptref(:,1),pierptref(:,2),lon0,lat0);

[pierptcomp(:,4),pierptcomp(:,5)] = absloc2relaloc(pierptcomp(:,1),pierptcomp(:,2),lon0,lat0);


diffdist = sqrt((pierptcomp(:,3)-pierptref(:,3)).^2 + (pierptcomp(:,4)-pierptref(:,4)).^2 + ...
    (pierptcomp(:,5)-pierptref(:,5)).^2);

if plotflag
    %% plot the map view
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
    
    scatter(f1.ax,tmr(:,6),tmr(:,7),2,[0.6 0.6 0.6],'filled','o');
    
    scatter(f1.ax,evt(:,8),evt(:,9),2,diffdist,'filled','o');
    % scatter(f1.ax,pierpt(:,1),pierpt(:,2),2,diffdist,'filled','o');
    
    % colormap(f1.ax,redblue);
    colormap(f1.ax,'jet');
    c=colorbar(f1.ax);
    c.Label.String = strcat('Abs. distance between piercing points (km)');
    c.Label.FontSize = 12;
    caxis(f1.ax,[0 2.5]);
    
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
    % patarea = [0 0;
    %            2.5 0;
    %            2.5 1;
    %            0 1;
    %            0 0];
    % patch(f5.ax,patarea(:,1),patarea(:,2),'k','Facealpha',0.1,'edgecolor','none');
    h = histogram(f2.ax,diffdist,'binwidth',0.5,'normalization','probability',...
        'facecolor','k','facea',0.5);
    perc95 = prctile(diffdist,95);
    
    for i = 1: length(h.Values)
        text(f2.ax,h.BinEdges(i)+0.1,h.Values(i),sprintf('%.2f',h.Values(i)));
    end
    text(f2.ax,0.5,0.95,num2str(size(diffdist,1)),'fontsize',13,'unit','normalized');
    text(f2.ax,0.3,0.85,sprintf('95 perc. %.2f',perc95),'fontsize',13,'unit','normalized');
    xlim(f2.ax,[0 5]);
    ylim(f2.ax,[0 1]);
    xlabel(f2.ax,'Abs. distance between piercing points (km)');
    ylabel(f2.ax,'Percentage');
    hold(f2.ax, 'off');
    
else
    f1=[];
    f2=[];
end

% keyboard


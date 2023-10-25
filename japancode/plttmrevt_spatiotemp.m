function [ax,evtobj,tmrobj] = plttmrevt_spatiotemp(regflag,evtall,tmrall,idbg,ided,magtol,dttol,dloctol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to plot the location of tremor detections and regular
% events during the same time period color-coded by time 
% 
%   INPUT:  
%       regflag: flag to indicate the region, 1 for western shikoku, 2 for kii
%                penninsula
%       slab:   slab grid file
%       evtall:  entire regular event catalog
%       tmrall:  entire tremor catalog
%       idbg:   begin day of the analysed period, counting from the start day
%       ided:   end day of the analysed period, counting from the start day
%       magtol:     magnitude threshold
%   OUTPUT:
%       ax:     current axes
%       evtets: event in that period
%       treets: tremor in that period
%       
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/18
% Last modified date:   2020/02/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % obtain the depth grid file for contour
% lat = 31:0.01:36;
% lon = 131:0.01:138;
% [longrd, latgrd] = meshgrid(lon,lat);
% depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd, 'linear');

d1 = datetime(tmrall(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
dateetsstr = datetime(d1+caldays(idbg:1:ided)-1,'Format','yyyy-MM-dd');
dateets = yyyymmdd(dateetsstr');

lend = length(dateets);
tmrets = [];
evtets = [];
for j = 1: lend
    tmrets = [tmrets; tmrall(dateets(j) == tmrall(:,1), :)];
    evtets = [evtets; evtall(dateets(j) == evtall(:,1) & evtall(:,11)>=magtol, :)];
end

di = datetime(tmrets(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
dtmrrela = caldays(between(d1,di,'days'))+1-idbg + tmrets(:,5)/24;

di = datetime(evtets(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
devtrela = caldays(between(d1,di,'days'))+1-idbg + evtets(:,5)/24;

%%% find the tremor and regular events that is close both in time and space
evtobj = [];
tmrobj = [];
for i = 1:size(evtets,1)
    tmrtmp = tmrets(abs(dtmrrela-devtrela(i))*24<=dttol, :);   % close in time
    if ~isempty(tmrtmp)
        dist = sqrt((tmrtmp(:,6)-evtets(i,8)).^2 + (tmrtmp(:,7)-evtets(i,9)).^2);
        tmrtmp2 = tmrtmp(dist<=dloctol, :);
        if ~isempty(tmrtmp2)
            evtobj = [evtobj; evtets(i,:)];
            tmrobj = [tmrobj; tmrtmp2];
        end
    end
end

if ~isempty(evtobj)
    disp(strcat(string(dateetsstr(1)),'to',string(dateetsstr(end))))
end

ax = gca;
hold(ax,'on');
if regflag == 1  % means western shikoku
    plot(ax,[134.1 134.7 132.11 131.51 134.1],[35.1 33.7 32.59 33.99 35.1],'color',[0.6 0.6 0.6],...
        'linew',2);
elseif regflag == 2  % means kii pen
    plot(ax,[134.6 135.2 137.0 136.4 134.6],[34.1 33.2 34.4 35.3 34.1],'color',[0.6 0.6 0.6],...
        'linew',2);
end
axis(ax, 'equal');
if regflag == 1
    xlim(ax,[131 135]);
    ylim(ax,[32.3 35.2]);
elseif regflag == 2
    xlim(ax,[134.5 137.2]);
    ylim(ax,[33 35.5]);
end
% contour(ax,longrd,latgrd,depgrd,10:10:100,'showtext','on','LabelSpacing',250,'linew',1,...
%     'linec',[0.8 0.8 0.8]);
coast = load('/home/data2/chaosong/matlab/Previous_mfiles/libBP/worldcoast.dat');
plot(ax,coast(:,1),coast(:,2),'black','linew',0.5);
p1=scatter(ax,tmrets(1,6),tmrets(1,7),6,[1 1 1],'filled','o','MarkerEdgeColor',...
    [0 0 0]);
p2=scatter(ax,evtets(1,8),evtets(1,9),6,[1 1 1],'filled','s','MarkerEdgeColor',...
    [0 0 0]);

scatter(ax,tmrets(:,6),tmrets(:,7),15,dtmrrela,'filled','o');
if ~isempty(tmrobj)
    scatter(ax,tmrobj(:,6),tmrobj(:,7),15,'o','MarkerEdgeColor',...
    [0 0 0]);
end

scatter(ax,evtets(:,8),evtets(:,9),40,devtrela,'filled','s');
if ~isempty(evtobj)
    scatter(ax,evtobj(:,8),evtobj(:,9),40,'s','MarkerEdgeColor',...
    [0 0 0]);    
end
colormap(ax,'jet');
c=colorbar(ax,'northoutside');
c.Label.String = strcat({'Days from '},string(dateetsstr(1)));
c.Label.FontSize = 10;
caxis(ax,[0 ided-idbg]);
legend(ax,[p1,p2],{ strcat('Tremors, n=',num2str(size(tmrets,1))), ...
    strcat('Regular >=',num2str(magtol),', n=',num2str(size(evtets,1))) },...
    'location','southeast');
ax.Box = 'on';
grid(ax, 'on');

hold(ax,'off');
% keyboard

%%% save figure
workpath = '/home/data2/chaosong/shikoku_kii';
figpath = strcat(workpath,'/figs');
if regflag == 1
    prefix = 'shikoku';
else
    prefix = 'kii';
end
if size(evtets,1) > 0
    print(gcf,'-dpdf',strcat(figpath,'/tmrevt.',prefix,'inets',...
          string(dateetsstr(1)),'to',string(dateetsstr(end)),'.pdf'));
end
% keyboard
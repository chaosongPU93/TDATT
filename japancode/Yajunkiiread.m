%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Yajunkiiread
%   read Yajun's the Kii Penninsula catalog and plot it, compare the depth 
%   information to the slab model in this area form USGS slab1.0, since
%   previously I realized that he might did something wrong with the slab
%   input for hypoinverse.
%
% NOTICE:
%   YAJUN's catalog is problematic, the time is reliable, but the lat,lon,dep
%   might be wrong to any extent, which can't be quantified now 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/12/22
% Last modified date:   2020/09/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear


set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not
scrsz=get(0,'ScreenSize');   % for example, scrsz=[1,1,1680,1050]
% scrsz=get(0,'MonitorPositions'); % return the position of all monitors

fid = fopen('/home/data2/chaosong/PAPERS/Peng_etal/Pengetal2018-sup-0003-supplementary.txt');
% fid = fopen('/home/data2/chaosong/PAPERS/Peng_etal/dummy');

%%%%%%%%%%%%% BIG NOTICE here!!! %%%%%%%%%%
% the original file copied from the web to linux actually has incompatible format issues at the end
% of each line with an extra '^M', which will cause malfunction of textscan, to remove all of them,
% in vim, type: :%s/\r//g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fspec = '%f %f %f %f';
datacell = textscan(fid, '%f %f %f %f \n','HeaderLines',1);
julday = floor(datacell{1}/3600/24)+1;
hour = floor((datacell{1}/3600/24+1 - julday)*24);
minute  = floor(((datacell{1}/3600/24+1 - julday)*24-hour)*60);
second = (((datacell{1}/3600/24+1 - julday)*24-hour)*60-minute)*60;
for i = 1: length(julday)
    if julday(i)<=366    % 2004 is leap
        year(i) = 2004;
        jday(i) = julday(i);
    elseif julday(i)<=366+365
        year(i) = 2005;
        jday(i) = julday(i)-366;
    elseif julday(i)<=366+365+365
        year(i) = 2006;
        jday(i) = julday(i)-366-365;
    elseif julday(i)<=366+365+365+365
        year(i) = 2007;
        jday(i) = julday(i)-366-365-365;
    elseif julday(i)<=366+365+365+365+366
        year(i) = 2008;
        jday(i) = julday(i)-366-365-365-365;
    elseif julday(i)<=366+365+365+365+366+365
        year(i) = 2009;
        jday(i) = julday(i)-366-365-365-365-366;
    elseif julday(i)<=366+365+365+365+366+365+365
        year(i) = 2010;
        jday(i) = julday(i)-366-365-365-365-366-365;
    elseif julday(i)<=366+365+365+365+366+365+365+365
        year(i) = 2011;
        jday(i) = julday(i)-366-365-365-365-366-365-365;
    elseif julday(i)<=366+365+365+365+366+365+365+365+366
        year(i) = 2012;
        jday(i) = julday(i)-366-365-365-365-366-365-365-365;
    elseif julday(i)<=366+365+365+365+366+365+365+365+366+365
        year(i) = 2013;
        jday(i) = julday(i)-366-365-365-365-366-365-365-365-366;
    elseif julday(i)<=366+365+365+365+366+365+365+365+366+365+365
        year(i) = 2014;
        jday(i) = julday(i)-366-365-365-365-366-365-365-365-366-365;
    elseif julday(i)<=366+365+365+365+366+365+365+365+366+365+365+365
        year(i) = 2015;
        jday(i) = julday(i)-366-365-365-365-366-365-365-365-366-365-365;
    elseif julday(i)<=366+365+365+365+366+365+365+365+366+365+365+365+366
        year(i) = 2016;
        jday(i) = julday(i)-366-365-365-365-366-365-365-365-366-365-365-365;
    end
    
    mo = jul2dat(year(i),jday(i));
    month(i) = mo(1);
    day(i) = mo(2);
end

rawcat = [datacell{1} year' jday' month' day' hour minute second datacell{2} datacell{3} datacell{4}];

%% plot the catalog, colorcoded by depth
f.fig = figure;
hold on
set(f.fig,'Position',[scrsz(1,1)*0.9 0 scrsz(1,3)*0.4 scrsz(1,4)*0.6]);
msize = 5;
scatter(rawcat(:,10), rawcat(:,9), msize, rawcat(:,11), 'filled','s');
colormap(jet)
colorbar
caxis([28, 33]);
axis([136.1 136.7 34.3 34.8]);
load coastlines
geoshow(coastlat,coastlon,'color','k');
plot([135 137 137 135 135], [33 33 35 35 33],'r-');
title("Yajun's 2018 Kii catalogue");
hold off

grd = load('/home/data2/chaosong/Seisbasics/hypoinverse/FromYajun/xyz.grid');
std = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/ryu_slab1.0_clip.xyz');   % most similar 


f.fig = figure;
set(f.fig,'Position',[scrsz(1,1)*0.9 0 scrsz(1,3)*0.4 scrsz(1,4)*0.6]);
scatter(grd(:,1),grd(:,2), 50, grd(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
% scatter(std(:,1),std(:,2), 50, -std(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([28, 33]);
axis([136.1 136.7 34.3 34.8]);
plot([135 137 137 135 135], [33 33 35 35 33],'r-');
geoshow(coastlat,coastlon,'color','k');
title("Corresponding slab depth at the same region");

% fid = fopen('/home/data2/chaosong/Seisbasics/hypoinverse/testhypo/slab1.0.grid','w+');
% fprintf(fid,'%.4f %.4f %.4f \n',ngrd');
% fclose(fid);








































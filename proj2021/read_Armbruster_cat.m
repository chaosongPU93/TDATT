 clc
clear
close all

workpath = getenv('ALLAN');

johncat = load(strcat(workpath,'/BOSTOCK/40sps.pgsssi.hyp'));


yr = floor(johncat(:,1)/10000);
yyyy = 2000+yr;
mm = floor((johncat(:,1)-yr*10000)/100);
dd = johncat(:,1)-yr*10000-mm*100;

hr = floor(johncat(:,2)/100);
mn = johncat(:,2)-hr*100;
time = hr*3600+mn*60+johncat(:,3);

lat = johncat(:,4)+johncat(:,5)/60;
lon = -(johncat(:,6)+johncat(:,7)/60);
dep = johncat(:,8);

[dx, dy] = absloc2relaloc(lon,lat,-123.585000, 48.436667);

newcat = [yyyy mm dd time dx dy lon lat dep];

dxloc0 = 0;
dyloc0 = 0;

%%
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

for i = 1: 12
  [x,y] = absloc2relaloc(str2num(stainfo(i,3)),str2num(stainfo(i,2)),...
    -123.585000, 48.436667);
  tmp = rad2deg(atan2(x, y));
  if tmp < 0
    staazi(i) = 360+tmp;
  else
    staazi(i) = tmp;
  end
end
disp(staazi([2 5 8 6]));  %klnb, pgc, ssib, silb
    
bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));    
nlfe = size(lfeloc,1);

%%% set path and name
catapath= ('/home/data2/chaosong/matlab/allan/BOSTOCK/');
cataname = ('total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(catapath, cataname));
%%% find the families that have LFE detections
famall = unique(catalog(:, 1));

%%% find the intersection between fam have detections and fam have locations
[~,~,ind] = intersect(famall, lfeloc(:, 1));
lfeuse = lfeloc(ind,:);
[dxlfe,dylfe] = absloc2relaloc(lfeuse(:,3),lfeuse(:,2),-123.585000, 48.436667);


%% plot
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
[scrsz, res] = pixelperinch(1);
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
nrow = 1;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    axis equal
end

% %%% reposition
% set(f2.ax(1), 'position', [ 0.08, 0.6, 0.36, 0.32]);
% set(f2.ax(2), 'position', [ 0.52, 0.6, 0.36, 0.32]);
% set(f2.ax(3), 'position', [ 0.08, 0.18, 0.36, 0.32]);
% set(f2.ax(4), 'position', [ 0.52, 0.18, 0.36, 0.32]);

xran = [-80 70];
yran = [-60 70];

msizehf = 4;
msizelf = 10;

% subplot 1 of figure i
hfplt = newcat;

ax = f.ax(1);
hold(ax,'on');
% scatter(ax,dxloc043,dyloc043,15,'k','filled');
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
[density1d,xloc2d,yloc2d,density2d] = density_matrix(hfplt(:,5),hfplt(:,6),...
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
% [density1d,xloc2d,yloc2d,density2d] = density_matrix(isolf(:,1),isolf(:,2),...
%     xran,yran,dxhf,dyhf);
% dum = density1d(density1d(:,3)>0, :);
% dum(dum(:,3)>1, :) = [];
% scatter(ax,dum(:,1),dum(:,2), 15, 'k','linew',1.5);
text(ax, 0.95, 0.93, "John's catalog",'FontSize',12,'unit','normalized','horizontalalignment','right',...
     'EdgeColor','k','Margin',2);
text(ax,0.85,0.8,'HF','FontSize',12,'unit','normalized','horizontalalignment','center');
text(ax,0.04,0.93,'a','FontSize',11,'unit','normalized','EdgeColor','k','Margin',2);
text(ax,0.4,0.93,strcat(num2str(length(hfplt)),{' detections'}),'FontSize',10,'unit','normalized',...
     'horizontalalignment','center');

%plot stations
for i=1:3
    % p1 is lzb trio stations
    p1 = plot(ax,dxlzb(i),dylzb(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','k');
    text(ax,dxlzb(i)+2,dylzb(i)+2,lzbtrio(i,1));
    % p2 is pgc trio stations
    p2 = plot(ax,dxpgc(i),dypgc(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','r');
    text(ax,dxpgc(i)+2,dypgc(i)+2,pgctrio(i,1)); 
end

for i = 1: length(otherind)
  % p3 is other unused stations
    p3 = plot(ax,dxother(i),dyother(i),'^','MarkerSize',8,'MarkerEdgeColor',...
        'k','markerf','w');
    text(ax,dxother(i)+2,dyother(i)+2,othersta(i,1)); 
end

scatter(ax,dxlfe, dylfe, 20, 'r', 'filled', 'MarkerEdgeColor', 'k');
text(ax,dxlfe, dylfe+1,num2str(lfeuse(:,1)), 'fontsize',12,'color','k');

ax.Box = 'on';
grid(ax, 'on');
axis(ax, 'equal');
ax.GridLineStyle = '--';
ax.XAxisLocation = 'top';
xlim(ax,xran);
ylim(ax,yran);
xticks(ax,xran(1):10:xran(2));
yticks(ax,yran(1):10:yran(2));
xlabel(ax,'E (km)','fontsize',11);
ylabel(ax,'N (km)','fontsize',11);
hold(ax,'off');

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forproj21');
print(f.fig,'-dpdf',strcat(getenv('MHOME'),'/Seisbasics/hypoinverse/forproj21',...
    '/johncatalogpgc.pdf'));


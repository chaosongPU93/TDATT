function compare_Yajun_Chao_Kii
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to compare the inverted locations of Yajun's and Chao's
% from the same time catalog from Yajun, as it seems that Yajun did something
% wrong with the slab grid file in hypoinverse 
% 
% 
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/09/20
% Last modified date:   2020/09/20
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% load grid file
close all
grd = load('/home/data2/chaosong/Seisbasics/hypoinverse/FromYajun/xyz.grid');
std = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/ryu_slab1.0_clip.xyz');   % most similar 

figure
scatter(grd(:,1),grd(:,2), 4, grd(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
plot([135 137 137 135 135], [33 33 35 35 33],'r-');
axis([135 139 32.5 35.5]);
title("Slab used by Yajun, color-coded by depth");

figure
scatter(std(:,1),std(:,2), 4, -std(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
plot([135 137 137 135 135], [33 33 35 35 33],'r-');
axis([135 139 32.5 35.5]);
title("Slab from USGS SLAB1.0, color-coded by depth");


ngrd = sortrows(grd, [2, 1], 'descend');

figure
subplot(3,1,1)
scatter(ngrd(:,1),ngrd(:,2), 4, ngrd(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
plot([136.1 136.7 136.7 136.1 136.1], [34.3 34.3 34.8 34.8 34.3],'r-');
axis([135 139 32.5 35.5]);
title("Yajun's slab, color-coded by depth");

fid = fopen('/home/data2/chaosong/Seisbasics/hypoinverse/FromYajun/new.xyz.grid','w+');
fprintf(fid,'%.4f %.4f %.4f \n',ngrd');
fclose(fid);


yajun = load('/home/data2/chaosong/Seisbasics/hypoinverse/FromYajun/eventloc.yajun');
chao = load('/home/data2/chaosong/Seisbasics/hypoinverse/FromYajun/eventloc.chao');


subplot(3,1,2)
scatter(-yajun(:,1),yajun(:,2), 4, yajun(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
% caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
% axis([135 139 32.5 35.5]);
caxis([28, 33]);
% axis([136.1 136.7 34.3 34.8]);
title("Yajun's locations, color-coded by depth");


subplot(3,1,3)
scatter(-chao(:,1),chao(:,2), 4, chao(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
% caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
% axis([135 139 32.5 35.5]);
caxis([28, 33]);
% axis([136.1 136.7 34.3 34.8]);
title("Chao's locations, color-coded by depth");


keyboard
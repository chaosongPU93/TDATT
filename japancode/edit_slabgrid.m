function edit_slabgrid
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to edit the slab1.0 grid file to make sure that it
% is compatible to hypoinverse
% 
% 
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/08/15
% Last modified date:   2019/08/15
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load grid file
% close all
grd = load('/home/data2/chaosong/Slab1.0_usgs/DepthGrid/cas_slab1.0_clip.xyz');
% std = load('/home/data2/chaosong/Seisbasics/hypoinverse/testhypo/old_xyz.grid');
std = load('/home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssitest/xyz.grid');


% lon = 121.92:0.01:125.76;
% nlon = length(lon);
% lat = 46.5:0.01:49.5;
% nlat = length(lat);
lon = 125.76:-0.01:121.92;
% lon = 125.5:-0.01:121.5;
nlon = length(lon);
lat = 49.5:-0.01:46.5;
nlat = length(lat);

[longrd, latgrd] = meshgrid(lon,lat);

lonmesh = reshape(longrd',[],1);
latmesh = reshape(latgrd',[],1);
% depnew = interp2(grd(:,1),grd(:,2),grd(:,3),lonmesh,latmesh,'spline');
F = scatteredInterpolant(grd(:,1),grd(:,2),grd(:,3),'natural','linear');
depmesh = F(lonmesh,latmesh);


%%% check the results
figure
subplot(3,1,1)
scatter(grd(:,1),grd(:,2), 4, grd(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
plot([125.76 125.76 121.92 121.92 125.76], [46.5 49.5 49.5 46.5 46.5],'r-');
axis([121 126 46 50]);
title('Slab 1.0, 0.02*0.02 deg');

subplot(3,1,2)
scatter(lonmesh, latmesh, 4, depmesh, 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
plot([125.76 125.76 121.92 121.92 125.76], [46.5 49.5 49.5 46.5 46.5],'r-');
axis([121 126 46 50]);
title('Slab 1.0, interpolated, 0.01*0.01 deg');

subplot(3,1,3)
scatter(std(:,1),std(:,2), 4, std(:,3), 'filled','o');   hold on  %, 'MarkerEdgeColor', 'w')
colormap(jet)
colorbar
caxis([round(min(grd(:,3))), round(max(grd(:,3)))]);
plot([125.76 125.76 121.92 121.92 125.76], [46.5 49.5 49.5 46.5 46.5],'r-');
axis([121 126 46 50]);
title('Old slab, 0.01*0.01 deg');

% %%% save file
% ngrd = [lonmesh latmesh depmesh];
% fid = fopen('/home/data2/chaosong/Seisbasics/hypoinverse/testhypo/slab1.0.grid','w+');
% fprintf(fid,'%.4f %.4f %.4f \n',ngrd');
% fclose(fid);



function rayinfo = evt_piercing_point(regflag,recalflag,slab,dslab,evtsel,stasel,kdense,velmod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to obtain the piercing point of the ray path from the 
% event to station on the slab interface. First, creat a kd-tree model 
% (KDtreeSearcher) based on the slab grid data (after interpolation). Second,
% use k-nearest neighbours (knn) searcher (knnsearch or a radius search 
% rangesearch) to find k-nearest points on the slab that are close to the 
% interpolated ray path. 
% 
%
% kdtree:
% 1.Kd-trees divide your data into nodes with at most BucketSize (default
%   is 50) points per node, based on coordinates (as opposed to categories).
%
%
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
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/09
% Last modified date:   2020/04/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
workpath = '/home/data2/chaosong/shikoku_kii';
datapath = strcat(workpath,'/matsave');

%% interpolate to get a much denser slab grid
% dslab = 0.002;
if regflag == 1
    prefix = 'shikoku';
    lat = 32.5:dslab:35.5;
    lon = 131.5:dslab:135;
elseif regflag == 2
    prefix = 'kii';
    lat = 32.5:dslab:35;
    lon = 134:dslab:137.5;
end
[longrd, latgrd] = meshgrid(lon,lat);

% create a interpolant object function for multiple usages
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');
depgrd = F(longrd,latgrd);
% depgrd = griddata(slab(:,1),slab(:,2),-slab(:,3),longrd,latgrd,'natural');

% % a much sparser one for plotting purpose only
% if regflag == 1
%     latspa = 32:0.02:35.5;
%     lonspa = 130.5:0.02:135.5;
% elseif regflag == 2
%     latspa = 32.5:0.02:35;
%     lonspa = 134:0.02:137.5;
% end
% [longrdspa, latgrdspa] = meshgrid(lonspa,latspa);
% depgrdspa = F(longrdspa,latgrdspa);
% % depgrdspa = griddata(slab(:,1),slab(:,2),-slab(:,3),longrdspa,latgrdspa,'natural');

% % for plotting purpose
% slabspa = [reshape(longrdspa,[],1) reshape(latgrdspa,[],1) reshape(depgrdspa,[],1)];
% [slabspac(:,1),slabspac(:,2)] = absloc2relaloc(slabspa(:,1),slabspa(:,2),lon0,lat0);
% slabspac(:,3) = slabspa(:,3);
% longrdspac = reshape(slabspac(:,1),size(longrdspa));
% latgrdspac = reshape(slabspac(:,2),size(latgrdspa));
% depgrdspac = reshape(slabspac(:,3),size(depgrdspa));


%% create a kdtree model of slab interface in cartesian coordinate
slabdense = [reshape(longrd,[],1) reshape(latgrd,[],1) reshape(depgrd,[],1)];
lat0 = 0.5*(lat(1)+lat(end));
lon0 = 0.5*(lon(1)+lon(end));
[slabcart(:,1),slabcart(:,2)] = absloc2relaloc(slabdense(:,1),slabdense(:,2),lon0,lat0);
slabcart(:,3) = slabdense(:,3);

kdtree = KDTreeSearcher(slabcart,'Distance','euclidean','BucketSize',50);


%% 
evt = evtsel;
sta = stasel;
nsta = size(sta,2);
nevt = size(evt,1);

% recalculation flag
% recalflag = 1;

% several available velocity models
if isequal(velmod,'jma2001')
    vmod = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.taup';
elseif isequal(velmod,'ukawa83')
    vmod = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ukawa83.taup';
elseif isequal(velmod,'ak135')
    vmod = 'ak135';
end

if recalflag
    % calculate the ray parameter between all stations and tremors
    tts = zeros(nevt, nsta);
    rayp = zeros(nevt, nsta);
    dist = zeros(nevt, nsta);
    pierpt = nan(nsta,nevt,3);
    edist = zeros(nevt, nsta);
    poorres = [];
    tic
    for i = 1: nsta
        stapt = zeros(nevt, 3);  % one event-sta pair--one piercing point with lon, lat, dep
        for j = 1: nevt
            % travel time and ray p
            tt=tauppath('mod',vmod,'dep',evt(j,10),'ph','s,S','evt',[evt(j,9) evt(j,8)],...
                        'sta',[sta(i).la sta(i).lo]);
            tts(j,i) = tt(1).time;   % choose the first S arrival
            rayp(j,i) = tt(1).rayparameter;  % unit in [s/deg]
            dist(j,i) = tt(1).distance;  % unit in [deg]
            raydep = tt(1).path.depth; % info of the ray path points, try 'help tauppath
            raylat = tt(1).path.latitude;
            raylon = tt(1).path.longitude;
            raypath = [raylon raylat raydep];
            
            raypathc = [];
            [raypathc(:,1),raypathc(:,2)] = absloc2relaloc(raypath(:,1),raypath(:,2),lon0,lat0);
            raypathc(:,3) = raypath(:,3);
            
%             if raypathc(end,2)>raypathc(1,2)
%                 dlat = 0.5;
%             else
%                 dlat = -0.5;
%             end
%             if raypathc(end,1)>raypathc(1,1)
%                 dlon = 0.5;
%             else
%                 dlon = -0.5;
%             end
%             nlat = floor((raypathc(end,2)-raypathc(1,2))/dlat)+1;
%             nlon = floor((raypathc(end,1)-raypathc(1,1))/dlon)+1;
%             npts = max(nlon,nlat);
%             londen = linspace(raypathc(1,1),raypathc(end,1),npts)';
%             latden = linspace(raypathc(1,2),raypathc(end,2),npts)';
            
%             kdense = 10;
            londen = kdensify1d(raypathc(:,1), kdense);
            latden = kdensify1d(raypathc(:,2), kdense);

            F2 = scatteredInterpolant(raypathc(:,1),raypathc(:,2),raypathc(:,3),'linear','linear');
            depden = F2(londen,latden);
            if ~isempty(depden)     % in case the triangulation is emty, points collinear
%             depden = griddata(raypathc(:,1),raypathc(:,2),raypathc(:,3),londen,latden,'linear');
                raydense = [londen latden depden];
                % eucdist is the min euclidean distance of each point in raydense to the kdtree
                % cluster
                [ind,eucdist] = knnsearch(kdtree,raydense,'K',1);
            else
                [ind,eucdist] = knnsearch(kdtree,raypathc,'K',1);
                poorres = [poorres; [i,j]];
            end
            
            % piercing point of that event on the slab interface
            % i.e. the one has the min eucdist, in fact the piercing point is an interpolated slab
            % point
            tmp = slabdense(ind(eucdist==min(eucdist)),:);
            stapt(j,:) = tmp;
            edist(j,i) = min(eucdist);

%             figure
%             ax=gca;
%             hold on;
%             box on;
%             grid on
%             surf(ax,longrdspac,latgrdspac,depgrdspac,'edgecolor','none');
%             plot3(ax,londen,latden,depden,'k-','linew',2);
%             plot3(ax,raypathc(:,1),raypathc(:,2),raypathc(:,3),'w--','linew',1.5);
%             [tmp(1),tmp(2)] = absloc2relaloc(tmp(1),tmp(2),lon0,lat0);
%             scatter3(ax,tmp(1),tmp(2),tmp(3),10,'r','filled','o');
%             set(ax,'ZDir','reverse');
%             colorbar;
%             view(ax,180,15);
%             xlabel('X (km)');
%             ylabel('Y (km)');
%             zlabel('Depth (km)');

        end
        
        pierpt(i,:,:) = stapt;
        
    end
    toc
    save(strcat(datapath,'/piercept_evt_bnd1_',num2str(dslab*100),'km_',num2str(kdense),'denseray_',...
         prefix,'_',velmod,'.mat'),'tts','rayp','dist','pierpt','poorres','edist');
    rayinfo.tts = tts ;
    rayinfo.rayp = rayp;
    rayinfo.dist = dist;
    rayinfo.pierpt = pierpt;
    rayinfo.poorres = poorres;
    rayinfo.edist = edist;
    
else
    rayinfo = load(strcat(datapath,'/piercept_evt_bnd1_',num2str(dslab*100),'km_',num2str(kdense),...
                'denseray_',prefix,'_',velmod,'.mat'));
%     rayinfo = load(strcat(datapath,'/piercept_evt_bnd2_',num2str(dslab*100),'km_',num2str(kdense),...
%                 'denseray_',prefix,'_',velmod,'.mat'));
end


% keyboard












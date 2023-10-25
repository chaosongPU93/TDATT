function color = colormatch(vect,cran,cmap,nbin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color = colormatch(vect,cmap,ncolor)
%
% In codes like 'scatter', you can easily use a data vector as the colorcoding
% by project the data point to the correct bin of the colormap. Here we want
% to do the similar thing, matching the each element in the 'vect' to the
% right bin of the 'cmap' that are separated into 'nbin' bins.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/12
% Last modified date:   2022/06/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('nbin',256);

if isequal(cmap, 'jet')
  cmat = jet(nbin);
elseif isequal(cmap, 'kelicol')
  load('kelim');  
  cmat = flipud(kelim);
  nbin = size(cmat,1);
elseif isequal(cmap, 'parula')
  cmat = parula(nbin);
elseif isequal(cmap, 'hsv')
  cmat = hsv(nbin);
elseif isequal(cmap, 'hot')
  cmat = hot(nbin);
elseif isequal(cmap, 'cool')
  cmat = cool(nbin);
elseif isequal(cmap, 'spring')
  cmat = spring(nbin);
elseif isequal(cmap, 'summer')
  cmat = summer(nbin);
elseif isequal(cmap, 'autumn')
  cmat = autumn(nbin);
elseif isequal(cmap, 'winter')
  cmat = winter(nbin);
elseif isequal(cmap, 'gray')
  cmat = gray(nbin);
elseif isequal(cmap, 'bone')
  cmat = bone(nbin);
elseif isequal(cmap, 'copper')
  cmat = copper(nbin);
elseif isequal(cmap, 'pink')
  cmat = pink(nbin);
elseif isequal(cmap, 'lines')
  cmat = lines(nbin);
elseif isequal(cmap, 'colorcube')
  cmat = colorcube(nbin);
elseif isequal(cmap, 'prism')
  cmat = prism(nbin);
elseif isequal(cmap, 'flag')
  cmat = flag(nbin);
elseif isequal(cmap, 'white')
  cmat = white(nbin);
end

binw = (cran(2)-cran(1))/nbin;

binedge = zeros(nbin,2);
for i = 1: nbin
  binedge(i,1) = cran(1)+ (i-1)*binw;
  binedge(i,2) = cran(1)+ i*binw;
  
end
bincnt = mean(binedge,2);

iwin = findwhichrange(vect,binedge);
iwin(iwin == -inf) = 1;
iwin(iwin == inf) = nbin;

color = cmat(iwin,:);

% keyboard

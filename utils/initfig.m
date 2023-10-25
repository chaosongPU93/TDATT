function [f] = initfig(widin,htin,nrow,ncol,ifig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = initfig(widin,htin,nrow,ncol,ifig)
%
% This function would create a figure composed by subplots with 'nrow' and 
% 'ncol'. Figure dimension is set by a width [widin] and height [htin] in 
% inches. Return the figure handle and axis handles as the attributes
% of 'f'.
%
% For A4 paper size, maximum width allowed is 8.5 inches, maximum height
% allowed is 11 inches. You definitely exceed this dimension, but printed
% figure to .pdf will be truncated to fit the paper, or you can resize it
% before printing.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/12
% Last modified date:   2022/03/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('widin',6);
defval('htin',5);
defval('nrow',1);
defval('ncol',1);
defval('ifig',[]);

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, resol] = pixelperinch(1);

if isempty(ifig)
  f.fig = figure;
else
  f.fig = figure(ifig);
end
f.fig.Renderer = 'painters';
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);

for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
    f.ax(isub).Box = 'on';
%     grid(f.ax(isub),'on');
end
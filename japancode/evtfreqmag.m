function f=pltevtmagfreq(f,evtmag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to analyse the relation between frequency and magnitude
% of the regular event catalog, spatial distribution or so
% 
%   INPUT:  
%       regflag: flag to indicate the region, 1 for western shikoku, 2 for kii
%                penninsula
%       depran: max depth range that allows the regular earthquakes fall below
%               the slab interface
%       recalflag: flag to indicate whether to re-cut the catalog
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/17
% Last modified date:   2020/02/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% plot the entire event catalog inside the study area, to find the suitable magnitude threshold
figure

histogram(f,evtmag,'binwidth',0.5,'normalization','count','facecolor','r','facea',0.8);
box on;
grid on;
mag75 = prctile(evtmag,75);
mag80 = prctile(evtmag,80);
mag85 = prctile(evtmag,85);
mag90 = prctile(evtmag,90);
text(f,0.6,0.9,strcat({'75 percentile: '},num2str(mag75)),'unit','normalized');
text(f,0.6,0.8,strcat({'80 percentile: '},num2str(mag80)),'unit','normalized');
text(f,0.6,0.7,strcat({'85 percentile: '},num2str(mag85)),'unit','normalized');
text(f,0.6,0.6,strcat({'90 percentile: '},num2str(mag90)),'unit','normalized');

return
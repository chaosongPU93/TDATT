function [ax,h] = pltevtmagfreq(ax,evtmag,binw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to analyse the relation between frequency and magnitude
% of the regular event catalog, spatial distribution or so
% 
%   INPUT:
%       evtmag:      magnitudes of all events
%       binw:     bin width
%
%   OUTPUT:
%       ax:    current axes
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/18
% Last modified date:   2020/02/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ax = gca;
hold(ax, 'on');
h = histogram(ax,evtmag','binwidth',binw,'normalization','count','facecolor','r','facea',0.8);
ax.Box = 'on';
grid(ax, 'on');
mag75 = prctile(evtmag,75);
mag80 = prctile(evtmag,80);
mag85 = prctile(evtmag,85);
mag90 = prctile(evtmag,90);
text(ax,0.65,0.9,strcat({'75 perc: '},sprintf('%.1f',mag75)),'unit','normalized');
text(ax,0.65,0.8,strcat({'80 perc: '},sprintf('%.1f',mag80)),'unit','normalized');
text(ax,0.65,0.7,strcat({'85 perc: '},sprintf('%.1f',mag85)),'unit','normalized');
text(ax,0.65,0.6,strcat({'90 perc: '},sprintf('%.1f',mag90)),'unit','normalized');
xlim(ax,[min(evtmag)-binw max(evtmag)+binw]);
xlabel(ax,'Magnitude');
ylabel(ax,'Counts');
hold(ax, 'off');

function f=pltevtduringtmr(f,d1,d2,ievtday,nevtobj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to plot the regular event counts against time during tremor
% active dates
%
%
%   INPUT:
%       f:      the construct containing frame and for figure
%       d1:     the start day of the time axis, Must be datestring format
%       d2:     the end day of the time axis, Must be datestring format
%       ievtday:    the regular event day on the time axis, counting from the start 
%       nevtobj:    the event counts on that day
%
%   OUTPUT:
%       f:    the construct for figure
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/17
% Last modified date:   2020/02/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndays = caldays(between(d1,d2,'days'))+1;

tkday = 1:ndays;
tkdaystr = string(datetime(d1+caldays(tkday)-1,'Format','yyyy-MM-dd'));

hold(f.ax(1),'on');
f.ax(1).Box = 'on';
grid(f.ax(1), 'on');
bar(f.ax(1),ievtday,nevtobj','Facec','r','Edgec','r');
f.ax(1).XTick = tkday(1:100:end);
f.ax(1).XTickLabel = tkdaystr(1:100:end);
f.ax(1).XTickLabelRotation = 45;
xlim(f.ax(1),[0 1500]);
ylabel(f.ax(1),'Counts');
title(f.ax(1),'Regular earthquakes during tremor active dates');
hold(f.ax(1),'off');


hold(f.ax(2),'on');
f.ax(2).Box = 'on';
grid(f.ax(2), 'on');
bar(f.ax(2),ievtday,nevtobj','Facec','r','Edgec','r');
f.ax(2).XTick = tkday(1:100:end);
f.ax(2).XTickLabel = tkdaystr(1:100:end);
f.ax(2).XTickLabelRotation = 45;
xlim(f.ax(2),[1501 3000]);
ylabel(f.ax(2),'Counts');
hold(f.ax(2),'off');

hold(f.ax(3),'on');
f.ax(3).Box = 'on';
grid(f.ax(3), 'on');
bar(f.ax(3),ievtday,nevtobj','Facec','r','Edgec','r');
f.ax(3).XTick = tkday(1:100:end);
f.ax(3).XTickLabel = tkdaystr(1:100:end);
f.ax(3).XTickLabelRotation = 45;
xlim(f.ax(3),[3001 4500]);
ylabel(f.ax(3),'Counts');
hold(f.ax(3),'off');

hold(f.ax(4),'on');
f.ax(4).Box = 'on';
grid(f.ax(4), 'on');
bar(f.ax(4),ievtday,nevtobj','Facec','r','Edgec','r');
f.ax(4).XTick = tkday(1:100:end);
f.ax(4).XTickLabel = tkdaystr(1:100:end);
f.ax(4).XTickLabelRotation = 45;
xlim(f.ax(4),[4501 6000]);
ylabel(f.ax(4),'Counts');
hold(f.ax(4),'off');





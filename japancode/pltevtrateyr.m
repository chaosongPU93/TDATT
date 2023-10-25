function f=pltevtrateyr(f,yrall,rtayr,rtiyr,rtoyr,rttmryr,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to plot the regular event rate per day per unit area inside
% or outside tremor regions, yearly variation, regardless of tremor occurring time
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/25
% Last modified date:   2020/02/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold(f.ax(1),'on');
f.ax(1).Box = 'on';
plot(f.ax(1),rtayr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtayr,1)
    plot(f.ax(1),rtayr(i,:),'o-','linew',1,'markers',4);
end
f.ax(1).XTick = 1:length(yrall);
f.ax(1).XTickLabel = num2str(yrall);
f.ax(1).XTickLabelRotation = 45;
xlabel(f.ax(1),'Year');
ylabel(f.ax(1),'Events per day per unit area');
xlim(f.ax(1),[1 length(yrall)]);
ylim(f.ax(1),[0 4]);
if class == 3
    ylim(f.ax(1),[0 7]);
end
title(f.ax(1),'Events in rectangle');
hold(f.ax(1),'off');

hold(f.ax(2),'on');
f.ax(2).Box = 'on';
plot(f.ax(2),rtiyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtiyr,1)
    plot(f.ax(2),rtiyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(2).XTick = 1:length(yrall);
f.ax(2).XTickLabel = num2str(yrall);
f.ax(2).XTickLabelRotation = 45;
xlabel(f.ax(2),'Year');
ylabel(f.ax(2),'Events per day per unit area');
xlim(f.ax(2),[1 length(yrall)]);
ylim(f.ax(2),[0 8]);
if class == 3
    ylim(f.ax(2),[0 7]);
end
title(f.ax(2),'Events inside tremor');
hold(f.ax(2),'off');

hold(f.ax(3),'on');
f.ax(3).Box = 'on';
plot(f.ax(3),rtoyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtoyr,1)
    plot(f.ax(3),rtoyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(3).XTick = 1:length(yrall);
f.ax(3).XTickLabel = num2str(yrall);
f.ax(3).XTickLabelRotation = 45;
xlabel(f.ax(3),'Year');
ylabel(f.ax(3),'Events per day per unit area');
xlim(f.ax(3),[1 length(yrall)]);
if class == 7
    ylim(f.ax(3),[0 4]);
    legend(f.ax(3),{'All mag','Mag < 0','Mag = 0-1','Mag = 1-2','Mag = 2-3','Mag = 3-4','Mag >= 4'}, ...
       'location','north' );
elseif class == 3
    ylim(f.ax(3),[0 7]);
    legend(f.ax(3),{'All mag','Mag < 1','Mag >= 1'},'location','north' );
end
title(f.ax(3),'Events outside tremor');
hold(f.ax(3),'off');   

hold(f.ax(4),'on');
f.ax(4).Box = 'on';
plot(f.ax(4),rttmryr,'ko-','linew',1,'markers',4);
f.ax(4).XTick = 1:length(yrall);
f.ax(4).XTickLabel = num2str(yrall);
f.ax(4).XTickLabelRotation = 45;
xlabel(f.ax(4),'Year');
ylabel(f.ax(4),'Tremors per day per unit area');
xlim(f.ax(4),[1 length(yrall)]);
ylim(f.ax(4),[4 8]);
legend(f.ax(4),{'Tremors'},'location','southwest');
title(f.ax(4),'Tremors inside tremor region');
hold(f.ax(4),'off');   
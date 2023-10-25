function f=pltevtratetmryr(f,yra,rtatyr,rtafyr,rtityr,rtifyr,rtotyr,rtofyr,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to plot the regular event rate per day per unit area during
% tremor or tremor free days, inside or outside tremor regions, yearly variation
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/25
% Last modified date:   2020/02/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold(f.ax(1),'on');
f.ax(1).Box = 'on';
plot(f.ax(1),rtatyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtatyr,1)
    plot(f.ax(1),rtatyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(1).XTick = 1:length(yra);
f.ax(1).XTickLabel = num2str(yra);
f.ax(1).XTickLabelRotation = 45;
xlabel(f.ax(1),'Year');
ylabel(f.ax(1),'Events per day per unit area');
xlim(f.ax(1),[1 length(yra)]);
ylim(f.ax(1),[0 6]);
if class == 3
    ylim(f.ax(1),[0 12]);
end
title(f.ax(1),'Events in rectangle in tremor days');
hold(f.ax(1),'off');


hold(f.ax(2),'on');
f.ax(2).Box = 'on';
plot(f.ax(2),rtafyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtafyr,1)
    plot(f.ax(2),rtafyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(2).XTick = 1:length(yra);
f.ax(2).XTickLabel = num2str(yra);
f.ax(2).XTickLabelRotation = 45;
xlabel(f.ax(2),'Year');
ylabel(f.ax(2),'Events per day per unit area');
xlim(f.ax(2),[1 length(yra)]);
ylim(f.ax(2),[0 6]);
if class == 3
    ylim(f.ax(2),[0 12]);
end
title(f.ax(2),'Events in rectangle in free days');
hold(f.ax(2),'off');


hold(f.ax(3),'on');
f.ax(3).Box = 'on';
plot(f.ax(3),rtityr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtityr,1)
    plot(f.ax(3),rtityr(i,:),'o-','linew',1,'markers',4);
end
f.ax(3).XTick = 1:length(yra);
f.ax(3).XTickLabel = num2str(yra);
f.ax(3).XTickLabelRotation = 45;
xlabel(f.ax(3),'Year');
ylabel(f.ax(3),'Events per day per unit area');
xlim(f.ax(3),[1 length(yra)]);
ylim(f.ax(3),[0 12]);
title(f.ax(3),'Events inside tremor in tremor days');
hold(f.ax(3),'off');


hold(f.ax(4),'on');
f.ax(4).Box = 'on';
plot(f.ax(4),rtifyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtifyr,1)
    plot(f.ax(4),rtifyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(4).XTick = 1:length(yra);
f.ax(4).XTickLabel = num2str(yra);
f.ax(4).XTickLabelRotation = 45;
xlabel(f.ax(4),'Year');
ylabel(f.ax(4),'Events per day per unit area');
xlim(f.ax(4),[1 length(yra)]);
ylim(f.ax(4),[0 12]);
title(f.ax(4),'Events inside tremor in free days');
hold(f.ax(4),'off');


hold(f.ax(5),'on');
f.ax(5).Box = 'on';
plot(f.ax(5),rtotyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtotyr,1)
    plot(f.ax(5),rtotyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(5).XTick = 1:length(yra);
f.ax(5).XTickLabel = num2str(yra);
f.ax(5).XTickLabelRotation = 45;
xlabel(f.ax(5),'Year');
ylabel(f.ax(5),'Events per day per unit area');
xlim(f.ax(5),[1 length(yra)]);
ylim(f.ax(5),[0 6]);
if class == 3
    ylim(f.ax(5),[0 12]);
end
title(f.ax(5),'Events outside tremor in tremor days');
hold(f.ax(5),'off');


hold(f.ax(6),'on');
f.ax(6).Box = 'on';
plot(f.ax(6),rtofyr(1,:),'ko-','linew',1,'markers',4);
for i=2:size(rtofyr,1)
    plot(f.ax(6),rtofyr(i,:),'o-','linew',1,'markers',4);
end
f.ax(6).XTick = 1:length(yra);
f.ax(6).XTickLabel = num2str(yra);
f.ax(6).XTickLabelRotation = 45;
xlabel(f.ax(6),'Year');
ylabel(f.ax(6),'Events per day per unit area');
xlim(f.ax(6),[1 length(yra)]);
if class == 7
    ylim(f.ax(6),[0 6]);
    legend(f.ax(6),{'All mag','Mag < 0','Mag = 0-1','Mag = 1-2','Mag = 2-3','Mag = 3-4','Mag >= 4'}, ...
        'location','north' );
elseif class == 3
    ylim(f.ax(6),[0 12]);
    legend(f.ax(6),{'All mag','Mag < 1','Mag >= 1'},'location','north' );
end
title(f.ax(6),'Events outside tremor in free days');
hold(f.ax(6),'off');











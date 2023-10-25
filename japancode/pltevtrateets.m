function f=pltevtrateets(f,rtaetsm,rtanom,rtietsm,rtinom,rtoetsm,rtonom,rttmrets,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to plot the regular event rate per day per unit area during
% ETS period or inter-ETS period, inside or outside tremor regions
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/26
% Last modified date:   2020/02/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold(f.ax(1),'on');
f.ax(1).Box = 'on';
plot(f.ax(1),rtaetsm(1,:),'k-','linew',1);
for i=2:size(rtaetsm,1)
    plot(f.ax(1),rtaetsm(i,:),'-','linew',1);
end
plot(f.ax(1),rttmrets/max(rttmrets)*5,'--','color',[0.5 0.5 0.5],'linew',1);
f.ax(1).XTick = 1:5:size(rtaetsm,2);
xlabel(f.ax(1),'ETS period number');
ylabel(f.ax(1),'Events per day per unit area');
xlim(f.ax(1),[1 size(rtaetsm,2)]);
ylim(f.ax(1),[0 6]);
if class == 3
    ylim(f.ax(1),[0 12]);
end
title(f.ax(1),'Events in rectangle during ETS');
hold(f.ax(1),'off');


hold(f.ax(2),'on');
f.ax(2).Box = 'on';
plot(f.ax(2),rtanom(1,:),'k-','linew',1);
for i=2:size(rtanom,1)
    plot(f.ax(2),rtanom(i,:),'-','linew',1);
end
f.ax(2).XTick = 1:5:size(rtanom,2);
xlabel(f.ax(2),'Inter ETS period number');
ylabel(f.ax(2),'Events per day per unit area');
xlim(f.ax(2),[1 size(rtanom,2)]);
ylim(f.ax(2),[0 6]);
if class == 3
    ylim(f.ax(2),[0 12]);
end
title(f.ax(2),'Events in rectangle during inter-ETS');
hold(f.ax(2),'off');


hold(f.ax(3),'on');
f.ax(3).Box = 'on';
plot(f.ax(3),rtietsm(1,:),'k-','linew',1);
for i=2:size(rtietsm,1)
    plot(f.ax(3),rtietsm(i,:),'-','linew',1);
end
plot(f.ax(3),rttmrets/max(rttmrets)*10,'--','color',[0.5 0.5 0.5],'linew',1);
f.ax(3).XTick = 1:5:size(rtietsm,2);
xlabel(f.ax(3),'ETS period number');
ylabel(f.ax(3),'Events per day per unit area');
xlim(f.ax(3),[1 size(rtietsm,2)]);
ylim(f.ax(3),[0 12]);
title(f.ax(3),'Events inside tremor during ETS');
hold(f.ax(3),'off');


hold(f.ax(4),'on');
f.ax(4).Box = 'on';
plot(f.ax(4),rtinom(1,:),'k-','linew',1);
for i=2:size(rtinom,1)
    plot(f.ax(4),rtinom(i,:),'-','linew',1);
end
f.ax(4).XTick = 1:5:size(rtinom,2);
xlabel(f.ax(4),'Inter ETS period number');
ylabel(f.ax(4),'Events per day per unit area');
xlim(f.ax(4),[1 size(rtinom,2)]);
ylim(f.ax(4),[0 12]);
title(f.ax(4),'Events inside tremor during inter-ETS');
hold(f.ax(4),'off');


hold(f.ax(5),'on');
f.ax(5).Box = 'on';
plot(f.ax(5),rtoetsm(1,:),'k-','linew',1);
for i=2:size(rtoetsm,1)
    plot(f.ax(5),rtoetsm(i,:),'-','linew',1);
end
plot(f.ax(5),rttmrets/max(rttmrets)*5,'--','color',[0.5 0.5 0.5],'linew',1);
f.ax(5).XTick = 1:5:size(rtoetsm,2);
xlabel(f.ax(5),'ETS period number');
ylabel(f.ax(5),'Events per day per unit area');
xlim(f.ax(5),[1 size(rtoetsm,2)]);
if class == 7
    ylim(f.ax(5),[0 6]);
    legend(f.ax(5),{'All mag','Mag < 0','Mag = 0-1','Mag = 1-2','Mag = 2-3','Mag = 3-4', ...
           'Mag >= 4','Norm. tremor'},'location','north' );
elseif class == 3
    ylim(f.ax(5),[0 12]);
    legend(f.ax(5),{'All mag','Mag < 1','Mag >= 1','Norm. tremor'},'location','north' );
end
title(f.ax(5),'Events outside tremor during ETS');
hold(f.ax(5),'off');


hold(f.ax(6),'on');
f.ax(6).Box = 'on';
plot(f.ax(6),rtonom(1,:),'k-','linew',1);
for i=2:size(rtonom,1)
    plot(f.ax(6),rtonom(i,:),'-','linew',1);
end
f.ax(6).XTick = 1:5:size(rtonom,2);
xlabel(f.ax(6),'Inter ETS period number');
ylabel(f.ax(6),'Events per day per unit area');
xlim(f.ax(6),[1 size(rtonom,2)]);
if class == 7
    ylim(f.ax(6),[0 6]);
    legend(f.ax(6),{'All mag','Mag < 0','Mag = 0-1','Mag = 1-2','Mag = 2-3','Mag = 3-4', ...
           'Mag >= 4'},'location','north' );
elseif class == 3
    ylim(f.ax(6),[0 12]);
    legend(f.ax(6),{'All mag','Mag < 1','Mag >= 1'},'location','north' );
end
title(f.ax(6),'Events outside tremor during inter-ETS');
hold(f.ax(6),'off');











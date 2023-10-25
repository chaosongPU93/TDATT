function varargout=axsym(ax,w,perc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = axsym(ax,w)
%
% This function is to make the axis range symmetric relative to 0, based on
% the maximum range of the current axis limit. This should be called after
% all data elements have been plotted. 
% --In case of any symbols being plotted right on the edge, here we expand
%   the range by 10 %. You can further expand the range by calling 'axranexp'.
% --If you know the proper axis range a priori, you should just specify it,
%   rather than using this automatic adjustment.
%
% INPUT:
% 
%  ax     Axis handles (default: gca)
%  w      1 for x-axes
%         2 for y-axes (default)
%         3 for x and y-axes
%  perc   expansion in percentage (default: 10)
%         
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/12
% Last modified date:   2022/03/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('ax',gca);
defval('w',2);
defval('perc',10);

for i = 1: length(ax)
  
  switch w
    case 1
      xran = (1+perc/100)*[-max(abs(ax(i).XLim)) max(abs(ax(i).XLim))];
      ax(i).XLim = xran;
    case 2
      yran = (1+perc/100)*[-max(abs(ax(i).YLim)) max(abs(ax(i).YLim))];
      ax(i).YLim = yran;
    case 3
      xran = (1+perc/100)*[-max(abs(ax(i).XLim)) max(abs(ax(i).XLim))];
      ax(i).XLim = xran;
      yran = (1+perc/100)*[-max(abs(ax(i).YLim)) max(abs(ax(i).YLim))];
      ax(i).YLim = yran;
  end
  
end

% Optional output
varns={ax};
varargout=varns(1:nargout);


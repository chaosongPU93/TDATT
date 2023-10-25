function windows = movingwins(indst,inded,mwlen,ovlplen,fplt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [windows] = movingwins(indst,inded,mwlen,ovlplen,fplt)
%
% This function is to automatically generate moving windows with a specified
% length 'mwlen' and a overlapping length 'ovlplen' between consecutive windows.
% The start and end indices for generation is given by 'indst' and 'inded'.
% The total number of moving windows is determined based on the inputs above.
% the function returns the start and end indices of each window, in the same
% context as 'indst' and 'inded'.
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/04/26
% Last modified date:   2022/04/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('fplt',0);   % default is not to plot resulting windows

%stepping length of moving wins
steplen = mwlen-ovlplen;

%total number of wins
nwin = floor((inded-indst+1-ovlplen)/steplen);

%compute the start and end indices for each window
windows = zeros(nwin, 2);
for i = 1: nwin
  istart = indst + (mwlen-ovlplen)*(i-1);
  iend = istart + mwlen -1;
  windows(i,1) = istart;
  windows(i,2) = iend;
end

%if some portions are not included
if windows(end,2) < inded
  %make it an separate window if its length is at least half of that of the others
  if inded-(windows(end,1)+(mwlen-ovlplen))+1 >= mwlen/2
    windows = [windows; [windows(end,1)+(mwlen-ovlplen) inded]];
    %otherwise just combine the portion to the last window
  else
    windows(end,2) = inded;
  end
end
      
%plot a schematic figure about all windows
if fplt
  figure
  hold on;
  box on; grid on
  for i = 1: nwin
    plot(windows(i,1): windows(i,2), i*ones(mwlen,1),'ko','MarkerFaceColor','b',...
      'MarkerEdgeColor','b','MarkerSize',4,'linew',1);
  end
  axis([0 inded+1 0 nwin+1]);
  plot([indst indst],[0 nwin+1], 'k--');
  plot([inded inded],[0 nwin+1], 'k--');
  text(0.1,0.9,sprintf('%d windows',nwin),'Units','normalized');
  xlabel('Index');
  ylabel('Window #');
end
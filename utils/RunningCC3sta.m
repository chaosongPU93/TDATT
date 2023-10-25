function [ircc,rcc12,rcc13,rcc23] = RunningCC3sta(trace,cclen,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is just to simplify the computation of RCC for 3-station pairs
% by calling 'RunningCC'. 'trace' needs to have 3 columns with the same order
% as in the defined station list
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/01/03
% Last modified date:   2023/01/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('method','cumsum');

[ircc,rcc12] = RunningCC(trace(:,1), trace(:,2), cclen, method);
[~,rcc13] = RunningCC(trace(:,1), trace(:,3), cclen, method);
[~,rcc23] = RunningCC(trace(:,2), trace(:,3), cclen, method);

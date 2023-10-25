%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the driving script to call 'mergeday_lf/hf_PGC' funcion to combine
% the detections of all days in any fam using PGC trio
%
%   The previous 'plotaddcheck_lf_PGC' is now served as plotting checking
%   results from addtional stations only because it will lose many 
%   information, and because in hypoinverse, the most useful result file
%   location is the detections with original offsets and time
%
% NOTE:
%   filtering effect correction should be
%   done in this step, i.e., set the fcflag to be 1
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/06/23
% Last modified date:   2021/06/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/PGCtrio');

nfampool = ['002';
            '243';
            '240';
            '253';
            '036';
            '251';
            ];        

        
%% for HF
for ifam = 1: size(nfampool,1) 
    fam = nfampool(ifam, :);
    fcflag = 1;  % set fcflag=1
    convflag = 0;   % convflag=0
    iup = 4;    % upsample number,times 
    mergeday_hf_PGC(fam,fcflag,convflag,iup);  

end


%% for LF
for ifam = 1: size(nfampool,1) 
    fam = nfampool(ifam, :);
    fcflag = 1;  % set fcflag=1
    convflag = 0;   % convflag=0
    iup = 1;    % upsample number,times
    mergeday_lf_PGC(fam,fcflag,convflag,iup);

end













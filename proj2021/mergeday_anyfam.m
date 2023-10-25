%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the driving script to call 'mergeday_lf/hf' funcion to combine
% the detections of all days in any fam using TWKB trio, similar to 
% 'mergeday_pgc_anyfam'
%
%
% NOTE:
%   filtering effect correction should be
%   done in this step, i.e., set the fcflag to be 1
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/24
% Last modified date:   2021/07/24
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
rstpath = strcat(datapath, '/LZBtrio');

% combine all fams into an pool, 12 fams in total, final version        
nfampool = ['002';
            '043';
            '141';
            '047';
            '010';
            '144';
            '099';
            '068';
            '125';
            '147';
            '017';
            '006';
            '001';
            ];
        
%% for HF
for ifam = 1: size(nfampool,1) 
    fam = nfampool(ifam, :);
    fcflag = 1;  % set fcflag=1
    convflag = 0;   % convflag=0
    iup = 4;    % upsample number,times 
    mergeday_hf(fam,fcflag,convflag,iup);  

end


%% for LF
for ifam = 1: size(nfampool,1) 
    fam = nfampool(ifam, :);
    fcflag = 1;  % set fcflag=1
    convflag = 0;   % convflag=0
    iup = 1;    % upsample number,times
    mergeday_lf(fam,fcflag,convflag,iup);

end













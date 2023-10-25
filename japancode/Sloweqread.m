function treall = Sloweqread(fpath, regflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to read the Obara NEID tremor catalog
% from papers: 1. Obara, K., Tanaka, S., Maeda, T., & Matsuzawa, T. (2010), GRL
% 2. Maeda, T., & Obara. K. (2009), JGR
%
% This function should also be capable to read other slow earthquake catalog
% from the Slow Earthquake Database http://www-solid.eps.s.u-tokyo.ac.jp/~sloweq/,
% because the data format is uniform. HAS NOT BEEN TESTED YET
%
%   NOTICE:
%       1. The converted catalog now using from also being found in NIED website
%           https://hinetwww11.bosai.go.jp/auth/tremor/auto_hypo_catalog/?_tm=201704&LANG=en
%       *2. The time zone is JST (UT+9), Japan Standard Time = UTC + 9 h. I have confirmed that
%           both catalogs are using the same time zone by looking at the same detection 
%
%   INPUT:  
%       fpath: the complete file path of the slow earthquakes/tremors from SED
%       regflag: flag to indicate the region, 1 for western shikoku, 2 for kii
%                penninsula
%   OUTPUT:
%       treall: the found slow earthquakes/tremors that meet the searching range 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/17
% Last modified date:   2020/02/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default value for easy debugging
defval('fpath', '/home/data2/chaosong/shikoku_kii/etscat/sloweq.20001231.6939.115617203.ObaraNEID.csv');
defval('regflag', 1);

% catalog data path
workpath = '/home/data2/chaosong/shikoku_kii';
tmrpath = strcat(workpath,'/etscat');

if regflag == 1
    prefix = 'shikoku';
elseif regflag == 2
    prefix = 'kii';
end


fid = fopen(fpath);
fspec = '%f %f %f %f %f %f %f %f %f %f %s %s %s \n';
datacell = textscan(fid, fspec, 'HeaderLines',23, 'Delimiter',',');
yr = datacell{1};
mo = datacell{2};
dy = datacell{3};
hr = datacell{4};
lat = datacell{7};
lon = datacell{8};

len = size(yr,1);
date = zeros(len,1);
for i = 1:len
    if mo(i) <= 9
        mostr = strcat('0',num2str(mo(i)));
    else
        mostr = num2str(mo(i));
    end
    
    if dy(i) <= 9
        dystr = strcat('0',num2str(dy(i)));
    else
        dystr = num2str(dy(i));
    end
    
    datestr = strcat(num2str(yr(i)),mostr,dystr);
    date(i,1) = str2num(datestr);
end
treall = [date yr mo dy hr lon lat];

% apply the region limit
if regflag == 1
    treall = treall(lat>=32.59 & lat<=35.1 & lon>=131.51 & lon<=134.7 & ...
                    3*lon-7*lat>=156.6 & 3*lon-7*lat<=168.2 & ...
                    7*lon+3*lat>=1022.54 & 7*lon+3*lat<=1044, :);
elseif regflag == 2  % means kii pen
    treall = treall(lat>=33.2 & lat<=35.3 & lon>=134.6 & lon<=137.0 & ...
                  2*lon-3*lat>=166.9 & 2*lon-3*lat<=170.8 & ...
                  3*lon+2*lat>=472 & 3*lon+2*lat<=479.8, :);
end

treall = sortrows(treall,[1, 5]);

% save to a matfile, in case of next use
save(strcat(tmrpath,'/tre',prefix,'.mat'),'treall');






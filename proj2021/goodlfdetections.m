function candidate = goodlfdetections(FALG,fam,nday)
% this stores the sequential number of qualified detections from 'way 2' in the 'identify_pgc.m'
% for each of the trial fams using the PGC trio, i.e., choosing the first 5 percent of detections
%
% 2003 062 is so glitchy!
% it seems like the total number of qualified isolated arrivals would be on the order of 10
%
% NOTE that the corresponding dates are:
% timoffrot= [
%             2003 061;
%             2003 062;
%             2003 063;
%             2004 196;
%             2004 197;
%             2004 198;
%             2004 199;
%             2005 254;
%             2005 255;
%             2005 256;
%             ];
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2021/07/01
% Last modified date:   2021/07/25

candidate = cell(1, nday);

if isequal(FALG,'PGC')
    if isequal(fam,'002')
        candidate{1} = [1 2];
        candidate{2} = [8];
        candidate{3} = [1 ];
        candidate{4} = [1]; %?
        candidate{5} = [1 2]; %?
        candidate{6} = [];
        candidate{7} = [1 2];
        candidate{8} = [];
        candidate{9} = [1 3 4 5 8];
        candidate{10} = [2 3 4 8];
        
    elseif isequal(fam,'243')
        candidate{1} = [];
        candidate{2} = [10];
        candidate{3} = [1 2 3 4];
        candidate{4} = [1]; %?
        candidate{5} = [1 3 10]; %?
        candidate{6} = [1 2 3];
        candidate{7} = [1 2 3 4 6 7 8 10];
        candidate{8} = [1 5];
        candidate{9} = [1 3 4 6 10 15];
        candidate{10} = [1 3 7 12 14 16 17];
        
    elseif isequal(fam,'253')
        candidate{1} = [1 4];
        candidate{2} = [];
        candidate{3} = [1];
        candidate{4} = [5]; %?
        candidate{5} = [1 7]; %?
        candidate{6} = [3];
        candidate{7} = [1 6];
        candidate{8} = [3 4 5];
        candidate{9} = [3];
        candidate{10} = [];
        
    elseif isequal(fam,'036')
        candidate{1} = [1 7];
        candidate{2} = [12];
        candidate{3} = [];
        candidate{4} = [3]; %?
        candidate{5} = [3]; %?
        candidate{6} = [];
        candidate{7} = [];
        candidate{8} = [];
        candidate{9} = [2];
        candidate{10} = [];
        
    elseif isequal(fam,'034')
        candidate{1} = [1];
        candidate{2} = [3 9];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = [7]; %?
        candidate{6} = [];
        candidate{7} = [1];
        candidate{8} = [4];
        candidate{9} = [1 2];
        candidate{10} = [];
        
    elseif isequal(fam,'061')
        candidate{1} = [2 3];
        candidate{2} = [1];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = []; %?
        candidate{6} = [];
        candidate{7} = [];
        candidate{8} = [];
        candidate{9} = [];
        candidate{10} = [];
        
    elseif isequal(fam,'023')
        candidate{1} = [1 3];
        candidate{2} = [];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = []; %?
        candidate{6} = [];
        candidate{7} = [];
        candidate{8} = [];
        candidate{9} = [];
        candidate{10} = [];
        
    elseif isequal(fam,'251')
        candidate{1} = [1 2 4];
        candidate{2} = [];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = []; %?
        candidate{6} = [];
        candidate{7} = [];
        candidate{8} = [];
        candidate{9} = [];
        candidate{10} = [];
        
    elseif isequal(fam,'240')
        candidate{1} = [3 4];
        candidate{2} = [];
        candidate{3} = [2 3];
        candidate{4} = [1 2 5]; %?
        candidate{5} = [1 5 6]; %?
        candidate{6} = [1 2];
        candidate{7} = [1 2];
        candidate{8} = [1 3];
        candidate{9} = [1 4 5 7 9];
        candidate{10} = [2 3 4 5];
        
    elseif isequal(fam,'255')
        candidate{1} = [];
        candidate{2} = [];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = [1]; %?
        candidate{6} = [];
        candidate{7} = [4];
        candidate{8} = [];
        candidate{9} = [];
        candidate{10} = [];
        
    elseif isequal(fam,'012')
        candidate{1} = [];
        candidate{2} = [];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = []; %?
        candidate{6} = [];
        candidate{7} = [];
        candidate{8} = [];
        candidate{9} = [];
        candidate{10} = [];
        
    elseif isequal(fam,'065')
        candidate{1} = [];
        candidate{2} = [];
        candidate{3} = [];
        candidate{4} = []; %?
        candidate{5} = []; %?
        candidate{6} = [5];
        candidate{7} = [1];
        candidate{8} = [];
        candidate{9} = [];
        candidate{10} = [1];
        
    end
    
elseif isequal(FALG,'TWKB')
    
    if isequal(fam,'147')
        %it seems like the total number of qualified isolated arrivals would be on the order of 10
        candidate = cell(1, nday);
        candidate{1} = [];
        candidate{2} = [];
        candidate{3} = [];
        candidate{4} = [1];
        candidate{5} = [];
        candidate{6} = [3]; %?
        candidate{7} = [];
        candidate{8} = [6];
        candidate{9} = [1 5 9 10 12 15 30]; %2
        candidate{10} = [6 8 21 26 30]; %?
        candidate{11} = [5 6 9 12 14 23 24 29 30]; %?
        candidate{12} = [2 8 11 17 ]; %15? 27?
        candidate{13} = []; %2? 5?
        candidate{14} = [4];
        candidate{15} = [15];   %?
        candidate{16} = [];
        candidate{17} = [6 7 ];
        candidate{18} = [5 7 8 10 11 14 19 20 22 25 26];    %37?
        candidate{19} = [14 24 26]; %3? 9? 11? 12? 13? 23?
        candidate{20} = [1 5 9];    %2? 3?
        candidate{21} = [8 ];   %2? 14? 15?
        candidate{22} = [1 43];   %15? 16? 18? 37? (1)39 47?
        candidate{23} = [1];    %4? 6? 10? 19? 26?
        
        
        
    end
    
end
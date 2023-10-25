function [evtintmrall,evtouttmrall] = Jmacatlogread2(regflag,regbound,depran,slab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to read the regular earthquake catalog from downloaded from
% JMA, https://www.data.jma.go.jp/svd/eqev/data/bulletin/hypo_e.html
%
% Version 2:
%   Compare to v1, now read in the boundary defined by tremor locations, and
%   return events inside or outside that boundary (now i am using the tremor active region
%   as the boundary, but any closed boundary is allowed)
%   
%   NOTICE:
%       1. Refer to the record format if not sure about anything
%           https://www.data.jma.go.jp/svd/eqev/data/bulletin/data/format/hypfmt_e.html
%       2. The time zone is JST (UT+9), Japan Standard Time = UTC + 9 h;
%   INPUT:  
%       regflag: flag to indicate the region, 1 for western shikoku, 2 for kii
%                penninsula
%       depran: max depth range that allows the regular earthquakes fall below
%               the slab interface
%       slab:   slab model
%   OUTPUT:
%       evtall: the found events that meet the searching range 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/24
% Last modified date:   2020/02/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % default value for easy debugging
% defval('regflag', 1);
% defval('depran', 10);

if regflag == 1
    prefix = 'shikoku';
elseif regflag == 2
    prefix = 'kii';
end

% catalog data path
workpath = '/home/data2/chaosong/shikoku_kii';
evtpath = strcat(workpath,'/jmacat02-18');


% time = 2001;
time1 = 2001:1:2017;
time = [time1'; 201801; 201802; 201803; 201804];

% define the interpolation object
F = scatteredInterpolant(slab(:,1),slab(:,2),-slab(:,3),'natural','linear');

evtintmrall = [];   % matrix
evtouttmrall = [];

eventintmrall = [];   % struct
eventouttmrall = [];
for it = 1: length(time)
    fid = fopen(strcat(evtpath,'/h',num2str(time(it))), 'r');
    format = '%96c \n';
    %     format = '%1c %4d %2d %2d %2d %2d %4.2f %4.2f %3d %4.2f %4.2f %4d %4.2f %4.2f %5.2f %3.2f %2.1f %1c %2.1f %1c %1c %1c %1c %1c %1c %1c %1d %3d %24c %3d %1c \n';
    %     format = '%1c %4d %2d %2d %2d %2d %4d %4d %3d %4d %4d %4d %4d %4d %5d %3d %2d %1c %2d %1c %1c %1c %1c %1c %1c %1c %1d %3d %24c %3d %1c \n';
    catcell = textscan(fid, format, [inf,1]);
    catstr = catcell{1};
    len = size(catstr,1);
    for i = 1: len
        if ~isempty(strtrim(catstr(i,1)))
            rtype(i,1) = catstr(i,1);   % record type
        else
            rtype(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,2:5)))
            yr(i,1) = str2num(catstr(i,2:5));   % year
        else
            yr(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,6:7)))
            mo(i,1) = str2num(catstr(i,6:7));   % month
        else
            mo(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,8:9)))
            dy(i,1) = str2num(catstr(i,8:9));   % day
        else
            dy(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,10:11)))
            hr(i,1) = str2num(catstr(i,10:11));     % hour
        else
            hr(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,12:13)))
            mi(i,1) = str2num(catstr(i,12:13));     % minute
        else
            mi(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,14:17)))
            sec(i,1) = str2num(catstr(i,14:17))/100;    % second
        else
            sec(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,18:21)))
            secse(i,1) = str2num(catstr(i,18:21))/100;      % standard error (se) of second
        else
            secse(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,22:24)))
            latdeg(i,1) = str2num(catstr(i,22:24));     % lat deg
        else
            latdeg(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,25:28)))
            latmin(i,1) = str2num(catstr(i,25:28))/100;     % lat min
        else
            latmin(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,29:32)))
            latminse(i,1) = str2num(catstr(i,29:32))/100;   % se of lat min
        else
            latminse(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,33:36)))
            londeg(i,1) = str2num(catstr(i,33:36));     % lon deg
        else
            londeg(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,37:40)))
            lonmin(i,1) = str2num(catstr(i,37:40))/100;     % lon min
        else
            lonmin(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,41:44)))
            lonminse(i,1) = str2num(catstr(i,41:44))/100;   % se of lon min 
        else
            lonminse(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,45:49)))
            if ~isempty(strtrim(catstr(i,48:49)))
                dep(i,1) = str2num(catstr(i,45:47))+str2num(catstr(i,48:49))/100;   % depth
            else
                dep(i,1) = str2num(catstr(i,45:49))/100;
            end
        else
            dep(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,50:52)))
            depse(i,1) = str2num(catstr(i,50:52))/100;  % se of depth
        else
            depse(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,53:54)))
            str = strtrim(catstr(i,53:54));
            if isequal(str(1),'-')
                mag1(i,1) = -0.1*str2num(str(2));   % magnitude 1
                
            elseif isequal(str(1),'A')
                mag1(i,1) = -1.0-0.1*str2num(str(2));
                
            elseif isequal(str(1),'B')
                mag1(i,1) = -2.0-0.1*str2num(str(2));
            elseif isequal(str(1),'C')
                mag1(i,1) = -3.0-0.1*str2num(str(2));
            else
                mag1(i,1) = str2num(str)/10;
            end
        else
            mag1(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,55)))
            mag1type(i,1) = catstr(i,55);   % magnitude 1 type
        else
            mag1type(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,56:57)))
            str = strtrim(catstr(i,56:57));
            if isequal(str(1),'-')
                mag2(i,1) = -0.1*str2num(str(2));   % magnitude 2
                
            elseif isequal(str(1),'A')
                mag2(i,1) = -1.0-0.1*str2num(str(2));
                
            elseif isequal(str(1),'B')
                mag2(i,1) = -2.0-0.1*str2num(str(2));
            elseif isequal(str(1),'C')
                mag2(i,1) = -3.0-0.1*str2num(str(2));
            else
                mag2(i,1) = str2num(str)/10;
            end
        else
            mag2(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,58)))
            mag2type(i,1) = catstr(i,58);   % magnitude 2 type
        else
            mag2type(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,59)))
            tttype(i,1) = catstr(i,59);     % travel time table type 
        else
            tttype(i,1) = 'N';      
        end
        
        if ~isempty(strtrim(catstr(i,60)))
            hyplocprec(i,1) = catstr(i,60);     % Hypocenter location precision
        else
            hyplocprec(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,61)))
            subinfo(i,1) = catstr(i,61);    % type of the event (1 for natural EQ; 5 for LFE)
        else
            subinfo(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,62)))
            maxinten(i,1) = catstr(i,62);   % maximum density
        else
            maxinten(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,63)))
            dmgcls(i,1) = catstr(i,63);     % damage class
        else
            dmgcls(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,64)))
            tsucls(i,1) = catstr(i,64);     % tsunami class
        else
            tsucls(i,1) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,65)))
            distrnum(i,1) = str2num(catstr(i,65));  % district number
        else
            distrnum(i,1) = 1e6;
        end
        
        
        if ~isempty(strtrim(catstr(i,66:68)))
            regnum(i,1) = str2num(catstr(i,66:68)); % region number
        else
            regnum(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,69:92)))
            regname(i,1:24) = catstr(i,69:92);      % region name
        else
            regname(i,1:24) = 'N';
        end
        
        if ~isempty(strtrim(catstr(i,93:95)))
            nsta(i,1) = str2num(catstr(i,93:95)); % num of sta for hypocenter determination
        else
            nsta(i,1) = 1e6;
        end
        
        if ~isempty(strtrim(catstr(i,96)))
            hypdetflag(i,1) = catstr(i,96);     % hypocenter quality flag
        else
            hypdetflag(i,1) = 'N';
        end
        
    end
    
    lat = latdeg+latmin/60;
    lon = londeg+lonmin/60;
    
    evtintmrper = [];   % matrix
    evtouttmrper = [];
    
    eventintmrper = [];   % struct
    eventouttmrper = [];
    % isinrect, logical value if the point is inside the predefined rectangle
    if regflag == 1  % means western shikoku
        % first a normal rectangle, then a tilted one indicate a belt
        isinrect = lat>=32.59 & lat<=35.1 & lon>=131.51 & lon<=134.7 & ...
            3*lon-7*lat>=156.6 & 3*lon-7*lat<=168.2 & ...
            7*lon+3*lat>=1022.54 & 7*lon+3*lat<=1044;
    elseif regflag == 2  % means kii pen
        isinrect = lat>=33.2 & lat<=35.3 & lon>=134.6 & lon<=137.0 & ...
            2*lon-3*lat>=166.9 & 2*lon-3*lat<=170.8 & ...
            3*lon+2*lat>=472 & 3*lon+2*lat<=479.8;
    end
    
    % isintmr, logical value if the point is inside the inner region defined by the boundary
    [is,ion] = inpolygon(lon,lat,regbound(:,1),regbound(:,2));
    isintmr = is | ion;
    
    for i = 1:len
                
        if isinrect(i)
            depslab = F(lon(i),lat(i));
            
            if dep(i)>=depslab && dep(i)<=depslab+depran
                
                evtyr = yr(i);
                evtmo = mo(i);
                evtdy = dy(i);
                evthr = hr(i);
                evtmin = mi(i);
                evtsec = sec(i);
                evtlat = lat(i);
                evtlon = lon(i);
                evtdep = dep(i);
                evtmag1 = mag1(i);
                evtmag2 = mag2(i);
                
                if evtmo <= 9
                    mostr = strcat('0',num2str(evtmo));
                else
                    mostr = num2str(evtmo);
                end
                
                if evtdy <= 9
                    dystr = strcat('0',num2str(evtdy));
                else
                    dystr = num2str(evtdy);
                end
                
                datestr = strcat(num2str(evtyr),mostr,dystr);
                evtdate = str2num(datestr);
                
                evt = [evtdate evtyr evtmo evtdy evthr evtmin evtsec evtlon evtlat ...
                    evtdep evtmag1 evtmag2];
                
                % a more complete struct array
                event.date = evtdate;
                event.yr = evtyr;
                event.mo = evtmo;
                event.dy = evtdy;
                event.hr = evthr;
                event.min = evtmin;
                event.sec = evtsec;
                event.secse = secse(i);
                event.lon = evtlon;
                event.lonminse = lonminse(i);
                event.lat = evtlat;
                event.latminse = latminse(i);
                event.dep = evtdep;
                event.depse = depse(i);
                event.mag1 = evtmag1;
                event.mag1tp = mag1type(i);
                event.mag2 = evtmag2;
                event.mag2tp = mag2type(i);
                event.tttype = tttype(i);
                event.hyplocprec = hyplocprec(i);
                event.evttp = subinfo(i);
                event.maxinten = maxinten(i);
                event.dmgcls = dmgcls(i);
                event.tsucls = tsucls(i);
                event.distrnum = distrnum(i);
                event.regnum = regnum(i);
                event.regnm = regname(i,:);
                event.nsta = nsta(i);
                event.hypdetflag = hypdetflag(i);
                
                if isintmr(i)
                    evtintmrper = [evtintmrper; evt];
                    eventintmrper = [eventintmrper; event];
                else
                    evtouttmrper = [evtouttmrper; evt];
                    eventouttmrper = [eventouttmrper; event];
                end
                                
            end
        end
    end
    evtintmrall = [evtintmrall; evtintmrper];
    evtouttmrall = [evtouttmrall; evtouttmrper];
    
    eventintmrall = [eventintmrall; eventintmrper];
    eventouttmrall = [eventouttmrall; eventouttmrper];
end

% save the matrix to a matfile, in case of next use
save(strcat(evtpath,'/regevt_in_tmrbound',prefix,'dep',num2str(depran),'.mat'),'evtintmrall');
save(strcat(evtpath,'/regevt_out_tmrbound',prefix,'dep',num2str(depran),'.mat'),'evtouttmrall');

% save the complete struct to a matfile, in case of next use
save(strcat(evtpath,'/regevtstruct_in_tmrbound',prefix,'dep',num2str(depran),'.mat'),...
     'eventintmrall');
save(strcat(evtpath,'/regevtstruct_out_tmrbound',prefix,'dep',num2str(depran),'.mat'),...
     'eventouttmrall');
















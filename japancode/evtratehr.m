function [rteirm,rteorm,rteotm,evtir,evtor,evtot,tmrir,tmror,tmrot] = ...
         evtratehr(evtinall,evtoutall,tmrall,areaall,areain,dttol,dloctol,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to calculate the event occurence rate per tremor detection
% per unit area that falls into the time range (dttol) relative to a tremor
% detection, inside a radius (dloctol), outside the radius, and outside the
% tremor boundary
%
% imagine there is a tremor detection at a time at a location, what is the detection rate inside
% a horizontal radius within a time range relative to the tremor timing; in the same time range,
% what is the rate for the region outside that circle, but inside tremor region boundary; And what
% is the case for the region outside tremor region boundary   
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/26
% Last modified date:   2020/02/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d1 = datetime(tmrall(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');

% relative time in days of all tremor detections
di = datetime(tmrall(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
dtmrrela = caldays(between(d1,di,'days')) + tmrall(:,5)/24;      % unit is day

% relative time in days of all inside tremor boundary events
di = datetime(evtinall(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
devtrela = caldays(between(d1,di,'days')) + (evtinall(:,5)+evtinall(:,6)/60+evtinall(:,6)/3600)/24;

% for the events inside or ouside the radius of one tremor detection, inside tremor boundary
evtir = [];
tmrir = [];
evtor = [];
tmror = [];

for i = 1:size(tmrall,1)
    evttmp = evtinall(abs(devtrela-dtmrrela(i))*24<=dttol, :);   % close in time
    if ~isempty(evttmp)
        dist = sqrt((evttmp(:,8)-tmrall(i,6)).^2 + (evttmp(:,9)-tmrall(i,7)).^2);
        evttmp2 = evttmp(dist<dloctol, :);
        evttmp3 = evttmp(dist>=dloctol, :);
        if ~isempty(evttmp2)
            evtir = [evtir; evttmp2];
            tmrir = [tmrir; tmrall(i,:)];
        end
        if ~isempty(evttmp3)
            evtor = [evtor; evttmp3];
            tmror = [tmror; tmrall(i,:)];
        end
    end
    
end

% events rate for inside the radius of one tremor detection and inside tremor boundary
evtir = unique(evtir,'rows','stable');
nevtinrad = size(evtir,1);
tmrir = unique(tmrir,'rows','stable');
ntmrinrad = size(tmrir,1);
[indr,arad] = boundary(tmrir(:,6),tmrir(:,7),0.5);
poly = polyshape(tmrir(indr,6),tmrir(indr,7));
peri = perimeter(poly);
arad = arad+dloctol*peri;
% narad = arad/areaall;
% narad = pi*dloctol^2/areaall;
narad = 1;

rtevtinrad = nevtinrad/ntmrinrad/narad;

if class == 7
    % events with different magnitudes
    nevtinradm0 = size(evtir(evtir(:,11)<0, :), 1);
    nevtinradm01 = size(evtir(evtir(:,11)>=0 & evtir(:,11)<1, :), 1);
    nevtinradm12 = size(evtir(evtir(:,11)>=1 & evtir(:,11)<2, :), 1);
    nevtinradm23 = size(evtir(evtir(:,11)>=2 & evtir(:,11)<3, :), 1);
    nevtinradm34 = size(evtir(evtir(:,11)>=3 & evtir(:,11)<4, :), 1);
    nevtinradm4 = size(evtir(evtir(:,11)>=4, :), 1);
    % occurence rate, num per tremor detection per unit area
    rtevtinradm0 = nevtinradm0/ntmrinrad/narad;
    rtevtinradm01 = nevtinradm01/ntmrinrad/narad;
    rtevtinradm12 = nevtinradm12/ntmrinrad/narad;
    rtevtinradm23 = nevtinradm23/ntmrinrad/narad;
    rtevtinradm34 = nevtinradm34/ntmrinrad/narad;
    rtevtinradm4 = nevtinradm4/ntmrinrad/narad;
    rteirm = [rtevtinrad; rtevtinradm0; rtevtinradm01; rtevtinradm12; rtevtinradm23; rtevtinradm34; ...
              rtevtinradm4];
    
elseif class == 3
    nevtinradml1 = size(evtir(evtir(:,11)<1, :), 1);
    nevtinradmg1 = size(evtir(evtir(:,11)>=1, :), 1);
    rtevtinradml1 = nevtinradml1/ntmrinrad/narad;
    rtevtinradmg1 = nevtinradmg1/ntmrinrad/narad;
    rteirm = [rtevtinrad; rtevtinradml1; rtevtinradmg1];
end


% events rate for otuside the radius of one tremor detection but inside tremor boundary
evtor = unique(evtor,'rows','stable');
nevtoutrad = size(evtor,1);
tmror = unique(tmror,'rows','stable');
ntmroutrad = size(tmror,1);
% naorad = (areain-pi*dloctol^2)/areaall;
naorad = 1;

rtevtoutrad = nevtoutrad/ntmroutrad/naorad;

if class == 7
    % events with different magnitudes
    nevtoutradm0 = size(evtor(evtor(:,11)<0, :), 1);
    nevtoutradm01 = size(evtor(evtor(:,11)>=0 & evtor(:,11)<1, :), 1);
    nevtoutradm12 = size(evtor(evtor(:,11)>=1 & evtor(:,11)<2, :), 1);
    nevtoutradm23 = size(evtor(evtor(:,11)>=2 & evtor(:,11)<3, :), 1);
    nevtoutradm34 = size(evtor(evtor(:,11)>=3 & evtor(:,11)<4, :), 1);
    nevtoutradm4 = size(evtor(evtor(:,11)>=4, :), 1);
    % occurence rate, num per tremor detection per unit area
    rtevtoutradm0 = nevtoutradm0/ntmroutrad/naorad;
    rtevtoutradm01 = nevtoutradm01/ntmroutrad/naorad;
    rtevtoutradm12 = nevtoutradm12/ntmroutrad/naorad;
    rtevtoutradm23 = nevtoutradm23/ntmroutrad/naorad;
    rtevtoutradm34 = nevtoutradm34/ntmroutrad/naorad;
    rtevtoutradm4 = nevtoutradm4/ntmroutrad/naorad;
    rteorm = [rtevtoutrad; rtevtoutradm0; rtevtoutradm01; rtevtoutradm12; rtevtoutradm23; ...
              rtevtoutradm34;rtevtoutradm4];

elseif class == 3
    nevtoutradml1 = size(evtor(evtor(:,11)<1, :), 1);
    nevtoutradmg1 = size(evtor(evtor(:,11)>=1, :), 1);
    rtevtoutradml1 = nevtoutradml1/ntmroutrad/naorad;
    rtevtoutradmg1 = nevtoutradmg1/ntmroutrad/naorad;
    rteorm = [rtevtoutrad; rtevtoutradml1; rtevtoutradmg1];
end
            
            
% relative time in days of all outside tremor boundary events
di = datetime(evtoutall(:,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
devtrela = caldays(between(d1,di,'days')) + (evtoutall(:,5)+evtoutall(:,6)/60+evtoutall(:,6)/3600)/24;

% for the events outside tremor boundary
evtot = [];
tmrot = [];

for i = 1:size(tmrall,1)
    evttmp = evtoutall(abs(devtrela-dtmrrela(i))*24<=dttol, :);   % close in time
    if ~isempty(evttmp)
        dist = sqrt((evttmp(:,8)-tmrall(i,6)).^2 + (evttmp(:,9)-tmrall(i,7)).^2);
        evttmp2 = evttmp(dist>=dloctol, :);
        if ~isempty(evttmp2)
            evtot = [evtot; evttmp2];
            tmrot = [tmrot; tmrall(i,:)];
        end
    end
end
            
% events rate for outside tremor boundary
evtot = unique(evtot,'rows','stable');
nevtouttmr = size(evtot,1);
tmrot = unique(tmrot,'rows','stable');
ntmrouttmr = size(tmrot,1);
% naotmr = (areaall-areain)/areaall;
naotmr = 1;

rtevtouttmr = nevtouttmr/ntmrouttmr/naotmr;

if class == 7
    % events with different magnitudes
    nevtouttmrm0 = size(evtot(evtot(:,11)<0, :), 1);
    nevtouttmrm01 = size(evtot(evtot(:,11)>=0 & evtot(:,11)<1, :), 1);
    nevtouttmrm12 = size(evtot(evtot(:,11)>=1 & evtot(:,11)<2, :), 1);
    nevtouttmrm23 = size(evtot(evtot(:,11)>=2 & evtot(:,11)<3, :), 1);
    nevtouttmrm34 = size(evtot(evtot(:,11)>=3 & evtot(:,11)<4, :), 1);
    nevtouttmrm4 = size(evtot(evtot(:,11)>=4, :), 1);
    % occurence rate, num per tremor detection per unit area
    rtevtouttmrm0 = nevtouttmrm0/ntmrouttmr/naotmr;
    rtevtouttmrm01 = nevtouttmrm01/ntmrouttmr/naotmr;
    rtevtouttmrm12 = nevtouttmrm12/ntmrouttmr/naotmr;
    rtevtouttmrm23 = nevtouttmrm23/ntmrouttmr/naotmr;
    rtevtouttmrm34 = nevtouttmrm34/ntmrouttmr/naotmr;
    rtevtouttmrm4 = nevtouttmrm4/ntmrouttmr/naotmr;
    rteotm = [rtevtouttmr; rtevtouttmrm0; rtevtouttmrm01; rtevtouttmrm12; rtevtouttmrm23; ...
              rtevtouttmrm34;rtevtouttmrm4];
          
elseif class == 3
    nevtouttmrml1 = size(evtot(evtot(:,11)<1, :), 1);
    nevtouttmrmg1 = size(evtot(evtot(:,11)>=1, :), 1);
    rtevtouttmrml1 = nevtouttmrml1/ntmrouttmr/naotmr;     
    rtevtouttmrmg1 = nevtouttmrmg1/ntmrouttmr/naotmr;      
    rteotm = [rtevtouttmr; rtevtouttmrml1; rtevtouttmrmg1];
end
          
          
          
          
          
          
          
          
          
          
          
          
          
          
            
 
            
            
            
            
            
            
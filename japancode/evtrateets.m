function [rtetsm,rtnom] = evtrateets(stday,edday,evtall,tmrall,narea,class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is function to calculate the event occurence rate per day per unit area
% in tremor days and tremor free days (depends on the dates fed) in every
% year within the input year range, and classify further according to the
% magnitude 
%   
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/25
% Last modified date:   2020/02/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d1  = datetime(tmrall(1,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
% d2  = datetime(tmrall(end,1),'InputFormat','yyyy-MM-dd','ConvertFrom','yyyymmdd');
% ndyall = caldays(between(d1,d2,'days'))+1;


rtnom = [];
for i = 1: length(stday)
    % start and end date of ETS 
    dis = datetime(d1+caldays(stday(i))-1,'Format','yyyy-MM-dd');
    die = datetime(d1+caldays(edday(i))-1,'Format','yyyy-MM-dd');
    dateis = yyyymmdd(dis');
    dateie = yyyymmdd(die');
    
    ndyets(i) = edday(i)-stday(i)+1;
    
    % all magnitude events
    evtets = evtall(evtall(:,1)>=dateis & evtall(:,1)<=dateie, :);
    nets(i) = size(evtets,1);
    rtets(i) = nets(i)/ndyets(i)/narea;

    if class == 7
        % events with different magnitudes
        netsm0(i) = size(evtets(evtets(:,11)<0, :), 1);
        netsm01(i) = size(evtets(evtets(:,11)>=0 & evtets(:,11)<1, :), 1);
        netsm12(i) = size(evtets(evtets(:,11)>=1 & evtets(:,11)<2, :), 1);
        netsm23(i) = size(evtets(evtets(:,11)>=2 & evtets(:,11)<3, :), 1);
        netsm34(i) = size(evtets(evtets(:,11)>=3 & evtets(:,11)<4, :), 1);
        netsm4(i) = size(evtets(evtets(:,11)>=4, :), 1);
        % occurence rate, num per day per unit area
        rtetsm0(i) = netsm0(i)/ndyets(i)/narea;
        rtetsm01(i) = netsm01(i)/ndyets(i)/narea;
        rtetsm12(i) = netsm12(i)/ndyets(i)/narea;
        rtetsm23(i) = netsm23(i)/ndyets(i)/narea;
        rtetsm34(i) = netsm34(i)/ndyets(i)/narea;
        rtetsm4(i) = netsm4(i)/ndyets(i)/narea;
        
    elseif class == 3
        netsml1(i) = size(evtets(evtets(:,11)<1, :), 1);
        netsmg1(i) = size(evtets(evtets(:,11)>=1, :), 1);
        rtetsml1(i) = netsml1(i)/ndyets(i)/narea;
        rtetsmg1(i) = netsmg1(i)/ndyets(i)/narea;
    end
    
    % events inter-ETS
    if i ~= length(stday)
        stno = edday(i);
        edno = stday(i+1);
    end
    dnos = datetime(d1+caldays(stno)-1,'Format','yyyy-MM-dd');
    dnoe = datetime(d1+caldays(edno)-1,'Format','yyyy-MM-dd');
    datenos = yyyymmdd(dnos');
    datenoe = yyyymmdd(dnoe');
    
    ndyno = edno-stno+1;
    % all magnitude events
    evtno = evtall(evtall(:,1)>=datenos & evtall(:,1)<=datenoe, :);
    nno = size(evtno,1);
    rtnoi = nno/ndyno/narea;
    
    if class == 7
        % events with different magnitudes
        nnom0 = size(evtno(evtno(:,11)<0, :), 1);
        nnom01 = size(evtno(evtno(:,11)>=0 & evtno(:,11)<1, :), 1);
        nnom12 = size(evtno(evtno(:,11)>=1 & evtno(:,11)<2, :), 1);
        nnom23 = size(evtno(evtno(:,11)>=2 & evtno(:,11)<3, :), 1);
        nnom34 = size(evtno(evtno(:,11)>=3 & evtno(:,11)<4, :), 1);
        nnom4 = size(evtno(evtno(:,11)>=4, :), 1);
        % occurence rate, num per day per unit area
        rtnom0i = nnom0/ndyno/narea;
        rtnom01i = nnom01/ndyno/narea;
        rtnom12i = nnom12/ndyno/narea;
        rtnom23i = nnom23/ndyno/narea;
        rtnom34i = nnom34/ndyno/narea;
        rtnom4i = nnom4/ndyno/narea;
        rtnomi = [rtnoi; rtnom0i; rtnom01i; rtnom12i; rtnom23i; rtnom34i; rtnom4i];
        rtnom = [rtnom rtnomi];
    
    elseif class == 3
        nnoml1 = size(evtno(evtno(:,11)<1, :), 1);
        nnomg1 = size(evtno(evtno(:,11)>=1, :), 1);
        rtnoml1i = nnoml1/ndyno/narea;
        rtnomg1i = nnomg1/ndyno/narea;
        rtnomi = [rtnoi; rtnoml1i; rtnomg1i];
        rtnom = [rtnom rtnomi];
    end
end


if class == 7 
    rtetsm = [rtets;rtetsm0;rtetsm01;rtetsm12;rtetsm23;rtetsm34;rtetsm4];
elseif class == 3
    rtetsm = [rtets;rtetsml1;rtetsmg1];
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
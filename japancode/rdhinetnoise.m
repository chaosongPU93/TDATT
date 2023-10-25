function [STAE,STAN,STAU,STAt] = ...
    rdhinetnoise(sacpath,stainfo,date,sps,lo,hi,npo,npa,pltflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to read the ambient nosie data from Japan Hinet
% (after or before instrument response removal), bandpass filtering, resample if
% needed. Return the whole day data
%
% 
%   NOTICE:
%       1. The time zone now used in both tremor, regular events catalog and the Hinet
%           data are JST (UT+9), Japan Standard Time = UTC + 9 h, which are consistent 
%           to each other. Unless specified, conversion to UTC frame is unnecessary
%
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/03/05
% Last modified date:   2020/03/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stnm = char(stainfo{1,1});    % station name, string

for i = 1: size(stnm,1)
% i=1;

    clear tmpe tmpn tmpu tmpef tmpnf tmpuf tmpefd tmpnfd tmpufd
    datestr = num2str(date);
    dirnm = datestr(3:end);
    staefnm = strcat(stnm(i,:),'.E');
    stanfnm = strcat(stnm(i,:),'.N');
    staufnm = strcat(stnm(i,:),'.U');
%     staefnm = strcat(stnm(i,:),'.E.new');
%     stanfnm = strcat(stnm(i,:),'.N.new');
%     staufnm = strcat(stnm(i,:),'.U.new');
    staepath = fullfile(sacpath,dirnm,staefnm);
    stanpath = fullfile(sacpath,dirnm,stanfnm);
    staupath = fullfile(sacpath,dirnm,staufnm);
    
    [tmpe,hdr,~,~,time]=readsac(staepath,0,'l');   %%% Function 'readsac'
    [tmpn,~,~,~,~]=readsac(stanpath,0,'l');
    [tmpu,~,~,~,~]=readsac(staupath,0,'l');
    
    spstrue = round(1./hdr.DELTA);

    %%% ALWAYS filter first, then resample, to avoid time-domain waveform distortion
    % filtering
    tmpef = Bandpass(tmpe,spstrue,lo,hi,npo,npa,'butter');
    tmpnf = Bandpass(tmpn,spstrue,lo,hi,npo,npa,'butter');
    tmpuf = Bandpass(tmpu,spstrue,lo,hi,npo,npa,'butter');
   
    % resample if needed
    [num, denom] = rat(sps/spstrue);
    tmpefd = resample(tmpef,num,denom);   % times = num/denom
    tmpnfd = resample(tmpnf,num,denom);
    tmpufd = resample(tmpuf,num,denom);
    time = linspace(round(time(1)), round(time(end)),...
                           round(length(time)*sps/spstrue))';
                           
    STAE(:,i) = tmpefd;     % 1-day E component data
    STAN(:,i) = tmpnfd;     % 1-day N component data
    STAU(:,i) = tmpufd;     % 1-day U component data
    STAt(:,i) = time;       % 1-day time data

end



if pltflag
    
    [scrsz, res] = pixelperinch(2);

    figure
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 11;   % maximum height allowed is 11 inches
    set(gcf,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    
    subplot(311)
    hold on
    for i = 1: size(stnm,1)
        plot(STAt(:,i),STAE(:,i)/max(STAE(:,i))+i,'color',[0.5 0.5 0.5]);
        text(size(STAE,1)*0.75,i,stnm(i,:),'fontsize',12);
    end
    title(strcat(datestr,{'  Normalized E component of selected ambient noise'}));    
    ylim([0 size(stnm,1)+1]);
    xlabel('Time (s)');
    box on
    hold off
    
    subplot(312)
    hold on
    for i = 1: size(stnm,1)
        plot(STAt(:,i),STAN(:,i)/max(STAN(:,i))+i,'color',[0.5 0.5 0.5]);
        text(size(STAN,1)*0.75,i,stnm(i,:),'fontsize',12);
    end
    title('Normalized N component of noise');
    ylim([0 size(stnm,1)+1]);
    xlabel('Time (s)');
    box on
    hold off
    
    subplot(313)
    hold on
    for i = 1: size(stnm,1)
        plot(STAt(:,i),STAU(:,i)/max(STAU(:,i))+i,'color',[0.5 0.5 0.5]);
        text(size(STAU,1)*0.75,i,stnm(i,:),'fontsize',12);
    end
    title('Normalized U component of noise');
    ylim([0 size(stnm,1)+1]);
    xlabel('Time (s)');
    box on
    hold off
    
end
function [evtbefhinet,evtnodata] = check_if_data_exist_func(evt,sacpath,datedown,evtbefhinet,...
                                                            evtnodata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to check if each query event that we are interested in
% is already in the hard drive;
% if yes, convert the data format from win32 to sac, remove the station response
% and so on if necessary;
% if not, output those dates and events to run the corresponding python script
% to download the data specific to each event;
% Also, if the dates are earlier than hinet online database, then still output
% those dates, but downloading those data needs to send online query to hinet
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/07/24
% Last modified date:   2020/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(evt(1),datedown)
    disp(evt(1))
    
    % output the select event date for the following bash and python scripts
    % note that this file contains only one date, so that the following operations are to deal
    % with one date one time, because each existing data trace has a length of one day, and it
    % takes considerable time to process
    fid = fopen(fullfile(sacpath,'selectdate'),'w+');
    fprintf(fid,'%d \n',evt(1));
    fclose(fid);
    
    % check if the raw sac data already exists, if not, convert the original hinet win32 data
    % to sac format first, i.e. filename like 'N.INMH.E'
    % based on single date and selected stations only to increase speed
    status = system('bash /home/data2/chaosong/japan_chao/datatosac_single.sh','-echo');
    
    % then check if the sac pz file for instrument response already exists, if not, try to
    % download those files as well
    % NOTE: it seems that python command is not compatible w/ matlab
    %          status = system('conda activate','-echo');
    %          status = system('python /home/data2/chaosong/japan_chao/extract_pz_single.py','-echo');
    
    % then check if the instrument response is already removed from the raw data,
    % i.e. filename like 'N.INMH.E.new',  if not, remove the station response
    status = system('bash /home/data2/chaosong/japan_chao/remove_response_single.sh','-echo');
    
else
    
    if evt(1) < 20040401    % hinet data starts from 2004-04-01
        fprintf('Target event date %d is earlier than Hinet database, needs to request online.\n',...
            evt(1));
        
        evtbefhinet = [evtbefhinet; evt];
        
    else
        fprintf('Target event date %d is not in hard drive, use down_hinet.py to request data.\n',...
            evt(1));
        
        evtnodata = [evtnodata; evt];
        
        % the following lines are commands for Python
        % this python script 'down_hinet.py' will download data & convert to sac & extract SAC
        % PZ file
        % NOTE: it seems that python command is not compatible w/ matlab
        %              status = system('conda activate','-echo');
        %              status = system('python /home/data2/chaosong/japan_chao/down_hinet.py','-echo');
        
        % the following lines are also not running here
        % this bash script 'remove_response_downhinet.sh' will remove the station response for
        % all data in the datelist
        %              status = system('bash /home/data2/chaosong/japan_chao/remove_response_downhinet.sh');
        
    end
end
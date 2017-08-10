function adcp_time_corr(adcpfile)
    
% function adcp_time_corr(adcpfile)
%
% Written to correct time lag of one hour caused by a wrong GPS
% setting (see logbook form July 2010 mission).
%
% F. Cyr, January 2012.
% ------------------------------------------------------- %
    
    
load(adcpfile)



if isfield(ADCP, {'timecorrected'}) == 1    
    disp('File already corrected, nothing to do!')
else
    
    ADCP.mtime = ADCP.mtime-1/24; % remove one hour
    ADCP.timecorrected = 1;
    save(adcpfile, 'ADCP')

end

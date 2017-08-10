%% Main header
command = sprintf('cat %s > %s', headerInfo, outfileODF);
system(command);

command = sprintf('sed -s "s/FILE_SPECIFICATION.*= ''''/FILE_SPECIFICATION = ''%s''/" %s > /tmp/tmp', outfileODF(1:end-4), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/EVENT_NUMBER.*= ''''/EVENT_NUMBER = ''%0.3d''/" %s > /tmp/tmp', eventCount, outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/CREATION_DATE.*= ''''/CREATION_DATE = ''%s''/g" %s > /tmp/tmp',datestr(now,0), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf(['sed -s "s/ORIG_CREATION_DATE.*= ''''/CREATION_DATE = ''%s''/g" %s > /tmp/tmp'],datestr(now,0), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf(['sed -s "s/START_DATE_TIME.*= ''''/START_DATE_TIME = ''%s''/" %s > /tmp/tmp'], datestr(mtime_eps(1),0), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/END_DATE_TIME.*= ''''/END_DATE_TIME = ''%s''/" %s > /tmp/tmp', datestr(mtime_eps(end),0), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf(['sed -s "s/INITIAL_LATITUDE.*= /INITIAL_LATITUDE = %3.5f/" %s > /tmp/tmp'], latVec(iprof), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/END_LATITUDE.*= /END_LATITUDE = %3.5f/" %s > /tmp/tmp', latVec(iprof), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf(['sed -s "s/INITIAL_LONGITUDE.*= /INITIAL_LONGITUDE = %3.5f/" %s > /tmp/tmp'], lonVec(iprof), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/END_LONGITUDE.*= /END_LONGITUDE = %3.5f/" %s > /tmp/tmp', lonVec(iprof), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/MIN_DEPTH.*= /MIN_DEPTH = %3.2f/" %s > /tmp/tmp', min(pres), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/MAX_DEPTH.*= /MAX_DEPTH = %3.2f/" %s > /tmp/tmp', max(pres), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf('sed -s "s/EVENT_NUMBER.*= /EVENT_NUMBER = ''%0.3d''/" %s > /tmp/tmp', iprof, outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

command = sprintf(['sed -s "s/DESCRIPTION.*= ''''/DESCRIPTION= ''%s''/" %s > /tmp/tmp'], archive_files, outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

%% Parameter headers:
% pressure

command = sprintf('cat %s >> %s', 'ODFheader_pres.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %3.2f,/" %s > /tmp/tmp', min(pres), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %3.2f,/" %s > /tmp/tmp', max(pres), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(pres));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(pres)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% pressure Flag
command = sprintf('cat %s >> %s', 'ODFheader_presFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %d,/" %s > /tmp/tmp', min(pres_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %d,/" %s > /tmp/tmp', max(pres_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(pres_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(pres_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% temperature
command = sprintf('cat %s >> %s', 'ODFheader_temp.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(temp), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(temp), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(temp));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(temp)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% temperature Flag
command = sprintf('cat %s >> %s', 'ODFheader_tempFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(temp_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(temp_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(temp_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(temp_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% conductivity
command = sprintf('cat %s >> %s', 'ODFheader_cond.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(cond), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(cond), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(cond));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(cond)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% conductivity Flag
command = sprintf('cat %s >> %s', 'ODFheader_condFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(cond_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(cond_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(cond_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(cond_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% salinity
command = sprintf('cat %s >> %s', 'ODFheader_sali.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(sali), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(sali), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(sali));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(sali)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% salinity Flag
command = sprintf('cat %s >> %s', 'ODFheader_saliFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(sali_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(sali_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(sali_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(sali_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% sigma-T
command = sprintf('cat %s >> %s', 'ODFheader_sigT.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(sigT), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(sigT), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(sigT));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(sigT)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% sigma-T Flag
command = sprintf('cat %s >> %s', 'ODFheader_sigTFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', min(sigT_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4f,/" %s > /tmp/tmp', max(sigT_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(sigT_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(sigT_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% Buoyancy freq.
command = sprintf('cat %s >> %s', 'ODFheader_buoy.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4e,/" %s > /tmp/tmp', min(buoy), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4e,/" %s > /tmp/tmp', max(buoy), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(buoy));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(buoy)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% Buoyancy freq. Flag
command = sprintf('cat %s >> %s', 'ODFheader_buoyFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %d,/" %s > /tmp/tmp', min(buoy_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %d,/" %s > /tmp/tmp', max(buoy_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(buoy_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(buoy_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% Dissipation rate of TKE.
command = sprintf('cat %s >> %s', 'ODFheader_diss.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4e,/" %s > /tmp/tmp', min(diss), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4e,/" %s > /tmp/tmp', max(diss), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(diss));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(diss)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% Dissipation rate of TKE Flag
command = sprintf('cat %s >> %s', 'ODFheader_dissFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %d,/" %s > /tmp/tmp', min(diss_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %d,/" %s > /tmp/tmp', max(diss_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(diss_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(diss_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% Diffusivity
command = sprintf('cat %s >> %s', 'ODFheader_diff.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %6.4e,/" %s > /tmp/tmp', min(diff), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %6.4e,/" %s > /tmp/tmp', max(diff), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(diff));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(diff)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

% Diffusivity Flag
command = sprintf('cat %s >> %s', 'ODFheader_diffFlag.txt', outfileODF);
system(command);
command = sprintf('sed -s "s/MINIMUM_VALUE= ,/MINIMUM_VALUE= %d/" %s > /tmp/tmp', min(diff_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/MAXIMUM_VALUE= ,/MAXIMUM_VALUE= %d,/" %s > /tmp/tmp', max(diff_flag), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);

I = find(isnan(diff_flag));
command = sprintf('sed -s "s/NUMBER_VALID= ,/NUMBER_VALID= %d,/" %s > /tmp/tmp', length(diff_flag)-length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);
command = sprintf('sed -s "s/NUMBER_NULL= ,/NUMBER_NULL= %d,/" %s > /tmp/tmp', length(I), outfileODF);
system(command);
command = sprintf('mv /tmp/tmp %s', outfileODF);
system(command);


%% Record header
command = sprintf('cat %s >> %s', 'ODFheader_record.txt', outfileODF);
system(command);







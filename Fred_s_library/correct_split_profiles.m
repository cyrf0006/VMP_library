% To be run in ~/WINDEX/data_processing/check_split_profiles
clear
repertories = 'slopeDays_july';

fid = fopen(repertories);
C = textscan(fid, '%s', 'delimiter', '\n');
rep = char(C{1});

nRep = size(rep, 1); %number of Repertories to look at


for i = 1:nRep
    
    repName = rep(i, :);
    I = find(repName==' ');   
    repName(I) = [];    
    disp(['Now in ' repName])
    
    eval(['mkdir ' repName(4:end)]);
    
    mypath = cd;
    eval(['cd ' repName(4:end)]);
    eval(['!cp /home/cyrf0006/WINDEX/data_processing/check_split_profiles/' repName '/DAT*.P ./']);
    eval(['!cp /home/cyrf0006/WINDEX/data_processing/check_split_profiles/' repName '/DAT.TXT ./']);
    eval(['!cp /home/cyrf0006/WINDEX/data_processing/check_split_profiles/' repName '/SETUP.TXT ./']);
    eval(['!cp /home/cyrf0006/WINDEX/data_processing/check_split_profiles/' repName '/CALIB.DAT ./']);

    
    % var_profile
    !ls -1 DAT*.P | sed 's/\.P//' > Pfiles
    disp('  -> var_profile')
    var_profile_cal('Pfiles')
    
    % eps_profile
    !ls -1 profile*.mat | sed 's/\.mat//' > prof_files
    disp('  -> eps_profile')
    eps_profile_fc('prof_files')

    % Back in orginal repertory
    eval(['cd ' mypath]);
end
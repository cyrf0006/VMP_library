function vmp_to_ODF(ctdFiles, epsFiles, gpsFiles, headerInfo)

% function vmp_to_ODF(ctdFiles, epsFiles, gpsFiles, headerInfo)
%
% This function will generate a series of *.ODF files (ex:
% VMP_YYYMMDD-hhmmss.ODF) 
%
%
% Inputs:
% 'ctdFiles' is a text file containing a list of path-to-files (ex:
% profilesXXX.mat). In linux, this is done for example with:
%      "ls -1 profile2011*.mat > windex_profiles2011.list"
% 'epsFiles' is a text file containing a list of path-to-files (ex: eps_profilesXXX.mat)      
%      "ls -1 eps_profile2011*.mat > windex_epsprofiles2011.list"
% 'gpsFiles' is a file containing GPS output in text. Here I merged
% all the gps files I had in a single one with shell command: 
%      $ cat 2012.txt >> mergedGPS
%      $ cat file2.txt >> mergedGPS,
% where mergedGPS should have the form:
%  48 28.728 68 30.677 2009 07 21 09 57 14 
%  48 28.728 68 30.677 2009 07 21 09 57 19 
%  48 28.728 68 30.676 2009 07 21 09 57 24 
%  (...)
%
% 'headerInfo' is a text file to-be-edited and containing something
% like:
%
%  ODF_HEADER,
%    FILE_SPECIFICATION = '',
%  CRUISE_HEADER,
%    COUNTRY_INSTITUTE_CODE = 18QO,
%    CRUISE_NUMBER = 'N/D',
%    ORGANIZATION = 'Ismer/MPO/Québec-Océan',
%    CHIEF_SCIENTIST = 'Frederic Cyr',
%    START_DATE = '21-JUN-2009 00:00:00.00',
%    END_DATE = '25-OCT-2012 23:59:00.00',
%    PLATFORM = 'Small craft boat (Macoma, Krill, BelugaI)',
%    CRUISE_NAME = 'WINDEX (2009-2012)',
%    CRUISE_DESCRIPTION = 'VMP casts in the Lower St. Lawrence Estuary (F. Cyr PhD
%   project)',
%  EVENT_HEADER,
%    DATA_TYPE= 'VMP',
%    EVENT_NUMBER= 'N/D',
%    EVENT_QUALIFIER1= '1',
%    EVENT_QUALIFIER2= 'DOWN',
%    CREATION_DATE= '',
%    START_DATE_TIME= '',
%    END_DATE_TIME= '',
%    INITIAL_LATITUDE= ,
%    INITIAL_LONGITUDE= ,
%    END_LATITUDE= ,
%    END_LONGITUDE= ,
%    MIN_DEPTH= ,
%    MAX_DEPTH= ,
%    SAMPLING_INTERVAL= N/D,
%    DEPTH_OFF_BOTTOM= N/D,
%    EVENT_COMMENTS1= 'Stations on section Rimouski',
%    EVENT_COMMENTS2= 'Data averaged in 1db vertical bins',
%  INSTRUMENT_HEADER,
%    INST_TYPE= 'Rockland Scientific - Vertical Microstructure profiler',
%    MODEL= 'VMP500',
%    SERIAL_NUMBER= '034 or 037',
%    INSTRUMENT_COMMENTS= 'EQUIPPED WITH SEA-BIRD CTD',
%  HISTORY_HEADER,
%    CREATION_DATE= '',
%    PROCESS= 'Conversion from .mat to ODF',
%
% 
% ** Note that this function also calls 'putODFheader.m', a script
% without input that must be manually edited to account for fields
% of the VMP files. This scripts also calls a series of text files
% (two file for each variable in the final ODF file), e.g. for the
% buoyancy frequency: ODFheader_buoy.txt and ODFheader_buoyFlag.txt
%
% usage ex: 
% in path-to-folder/toODF:
%   vmp_to_ODF('windex_profile2011.list', 'windex_epsprofile2011.list','mergedGPS','ODFheader2011.txt')


% Fred's personal examples: 
%in /media/Seagate Backup Plus Drive/BackupCarbon/WINDEX/data_processing/toODF:
%   vmp_to_ODF('tadoussac_profile.list', 'tadoussac_epsprofile.list','mergedGPS_tadoussac', 'ODFheader_tadoussac.txt')
%   vmp_to_ODF('windex_20110922_profile.list', 'windex_20110922_epsprofile.list','mergedGPS','ODFheader.txt')
%   vmp_to_ODF('windex_profile.list', 'windex_epsprofile.list','mergedGPS','ODFheader.txt')

% author: F. Cyr - june 2014
%   modified - Feb. 2015
%   modified - May 2015
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% few params
N2_min = -6;
N2_count = 0;
GAMMA = 0.2;



%% Get lat-lon of profiles
[latVec, lonVec] = whereare(epsFiles, gpsFiles);


% load *.P files names (file in which are recorded .P files)
fid = fopen(ctdFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
fineScale = char(C{1});

fid = fopen(epsFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
microScale = char(C{1});

noProfiles = size(fineScale, 1);
for iprof = 1:noProfiles
    
    %% Read files
    fname = fineScale(iprof, :);
    I = find(fname==' ');   
    fname(I) = [];
    load(fname);
    
    fname = microScale(iprof, :);
    I = find(fname==' ');   
    fname(I) = [];
    load(fname);
    
    %% Isolate archive names (DAT*.P)
    if exist('Pfile_archive')
        MATfile_archive = [Pfile_archive(1:end-1) 'mat'];
        archive_files = [Pfile_archive ' ' MARfile_archive];
    else % F. Cyr PhD data (2009-2013) will fall into this
        archive_files = ['DAT_' datestr(mtime_eps(1),29) '_0*.P' ' ' 'DAT_' datestr(mtime_eps(1),29) '_0*.mat'];        
    end
    
    
    %% Calculation of turbulence variables
    if length(p_eps1) == length(p_eps2) % no problem, easy
        p_k = p_eps1;
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        EPS2(:,2) = eps2;
    elseif min(p_eps1)<min(p_eps2) % p1 drive, nomatter size of p2!
        p_k = p_eps1;        
        p_N = p_eps1;
        EPS2=nan(length(p_eps1),2);
        EPS2(:,1) = eps1;
        ind = find(p_eps2 == p_eps1(1));
        EPS2(ind:length(p_eps2),2) = eps2;
    else % p2 drive, nomatter size of p1!
        p_k = p_eps2;
        p_N = p_eps1;
        EPS2=nan(length(p_eps2),2);
        EPS2(:,2) = eps2;
        ind = find(p_eps1 == p_eps2(1));
        EPS2(ind:length(p_eps1),1) = eps1;
    end
   
    % "selection" average
    if size(p_k, 1) < size(p_k, 2);
        p_k = p_k';
    end
    I = find(EPS2(:,1)>10*EPS2(:,2));
    EPS2(I,1)=NaN;
    I = find(EPS2(:,2)>10*EPS2(:,1));
    EPS2(I,2)=NaN;
   
    EPS = nanmean(EPS2,2); %MEAN 
    
    % Homestyle despike (only if not full of nans)
    if sum(~isnan(EPS)) ~= 0        
        [Y, No] = Dspike(EPS, 5, 8);
        
        % uses home despike
        EPS = Y;
    else
        disp('next cast will have only NaNs for turbulence:')
    end
        
    % mean N2
    N2=nan(length(p_k),1);
    if length(p_N)==length(p_k);
        N2(:) = N.^2;
    else
        N2(ind:length(p_N))=N.^2;
    end
    
    % compute K_rho
    K_rho = GAMMA.*EPS./(N2);

    % Remove unrealistic diffusivity
    I = find(log10(N2)<N2_min);
    if ~isempty(I)==1
        K_rho_ignored(N2_count+1:N2_count+length(I))=K_rho(I)';
        K_rho(I)=NaN;
        N2_count = N2_count+length(I);
    end
    
    %% Bin fine scale to turbulence resolution
    Tbin = nan(size(p_k));
    Cbin = nan(size(p_k));
    Sbin = nan(size(p_k));
    zbin = p_k(2)-p_k(1);
    for i = 1:length(p_k)
        I = find(P>=p_k(i)-zbin/2 & P<p_k(i)+zbin/2);
        Tbin(i) = nanmean(SBT(I));
        Cbin(i) = nanmean(SBC(I));
        Sbin(i) = nanmean(SBS(I));        
    end
    
     % sigma-t   
    SA = gsw_SA_from_SP(Sbin, p_k, -70, 48);
    CT = gsw_CT_from_t(Sbin, Tbin, p_k);
    sigT = gsw_sigma0(SA,CT);    
      
    %% Write output in an ascii file
    thedate = datestr(mtime(1),30);
    thedate(9) = '-';
    % outfileODF = ['CTD_' thedate '.ODF']; % previous outname

    if iprof == 1
        eventCount = 1;
        cruiseDate = thedate(1:8);
    elseif strcmp(cruiseDate, thedate(1:8))~=1
        cruiseDate = thedate(1:8);
        eventCount = 1;
    end
        
    %% Build filename
% $$$     command = sprintf('grep CRUISE_NUMBER %s | sed -s "s/\\(^.*CRUISE_NUMBER.*''\\)\\([0-9]\\{7\\}\\)\\(''.*\\,$\\)/\\2/"' ,headerInfo);
% $$$     [status,cruiseNo] = system(command);
% $$$     
% $$$     command = sprintf('grep EVENT_QUALIFIER1 %s | sed -s "s/\\(^.*EVENT_QUALIFIER1.*''\\)\\([0-9]\\{1\\}\\)\\(''.*,$\\)/\\2/"' ,headerInfo);
% $$$     [status,qual1] = system(command);
% $$$     
% $$$     command = sprintf('grep EVENT_QUALIFIER2 %s | sed -s "s/\\(^.*EVENT_QUALIFIER2.*''\\)\\([A-Z]\\{2\\}\\)\\(''.*,$\\)/\\2/"' ,headerInfo);
% $$$     [status,qual2] = system(command);
% $$$     
% $$$     cruiseNo(end) = []; % Remove line break (this I don't know why!)
% $$$     qual1(end) = [];
% $$$     qual2(end) = [];   
        
    % Here it is weak, but much simpler than the preceeding:
    qual1 = 1; %unless if file already exist (see below)
    qual2 = 'DN'; % always 'downcast' for VMP    
    
    outfileODF = sprintf('CTD_%s_%0.3d_%d_%s.ODF', cruiseDate, eventCount, qual1, qual2);    

    % rename fields
    pres = p_k;
    temp = Tbin;
    cond = Cbin;
    sali = Sbin;
    sigT = sigT;
    diss = EPS;
    diff = K_rho;
    buoy = N2;
    
    % flags
    pres_flag = zeros(size(pres));
    temp_flag = zeros(size(pres));
    cond_flag = zeros(size(pres));
    sali_flag = zeros(size(pres));
    sigT_flag = zeros(size(pres));
    diss_flag = zeros(size(pres));
    diff_flag = zeros(size(pres));
    buoy_flag = zeros(size(pres));

    I = find(isnan(temp));
    temp(I) = -99.0000;
    temp_flag(I) = 9;
    I = find(isnan(cond));
    cond(I) = -99.0000;
    cond_flag(I) = 9;
    I = find(isnan(sali));
    sali(I) = -99.0000;
    sali_flag(I) = 9;
    I = find(isnan(sigT));
    sigT(I) = -99.0000;
    sigT_flag(I) = 9;
    I = find(isnan(diss));
    diss(I) = -99.0000;
    diss_flag(I) = 9;
    I = find(isnan(diff));
    diff(I) = -99.0000;
    diff_flag(I) = 9;
    I = find(isnan(buoy));
    buoy(I) = -99.0000;
    buoy_flag(I) = 9;
    
    %% Build the ODF file
    disp(['saving ' outfileODF]);
    fid = fopen(outfileODF);
    if fid ~= -1
        outfileODF(end-7) = '2'; 
        disp([' File exist, create file ' outfileODF ' instead']);
        fclose(fid);
    end

    % output matrix in temporary file
    outputMatrix = [pres pres_flag temp temp_flag cond cond_flag ...
                    sali sali_flag sigT sigT_flag buoy buoy_flag diss diss_flag diff diff_flag];
    
    double(outputMatrix);
    fid = fopen('/tmp/VMP.dat', 'wt'); %temporary file 
    fprintf(fid, ['%10.4f %d %10.4f %d %10.4f %d %10.4f %d %10.4f %d %10.4e %d %10.4e %d %10.4e %d\n'], outputMatrix');
    fclose(fid);    
    
    % Put header
    putODFheader;
    
    % add data    
    command = sprintf('cat /tmp/VMP.dat >> %s', outfileODF);
    system(command);
    system('rm /tmp/VMP.dat');
    
    eventCount = eventCount+1;
    
end 
disp(sprintf('%d files created', noProfiles))
toc
%%%%%%%%%%%%%%%%%%%%
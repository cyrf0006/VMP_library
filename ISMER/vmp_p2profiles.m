function vmp_p2profiles(file_names, varargin);

% function vmp_p2profiles(file_names, varargin);
%
% Before 2013-07-26 was known as "var_profile_cal.m". Help menu is
% thus a copy-paste of the preceeding function.
%
%
% This script splits each variables of the day on separated profiles.
%
% The parameter is a file containing .P files that we want to split into
% profiles. For instance, if we have DAT000.P, DAT001.P and DAT002.P,
% file_names would be: file_names = [ DAT000
%                                     DAT001
%                                     DAT002 ]
% note that the extension must be removed.
% In linux, an easy command to do in folder containing *.P is:
%
% "ls -1 DAT*.P | sed 's/\.P//' > Pfiles.list"
%
% In this example, Pfiles would be the parameter of the function
% and the calling would be : vmp_p2profiles('Pfiles').
%
% It is worth to note that in addition of Pfiles, 2 other files
% must be in the current folder: SETUP.TXT and CALIB.DAT (see
% p2mat.m help menu for details...
%
% The function returns nothing, but first creates a .mat file for each .P
% file (using p2mat.m), and then creates a .mat file for each profile
% contained in those files. If the total number of profiles in files DAT000
% to DAT002 is 15, var_profile(Pfiles) will create "profile001.mat" to "profile015.mat"
%
% This script also calibrates microstructure salinity and temperature
% before saving them. A version of this script with no calibration is
% called var_profile.m but is not updated
%
% Finally, this script creates .dat files with the 4 fields needed by the
% Reorange method from Galbraith & Kelley 1996. The created files named
% fineXXX.dat and microXXX.dat for the fine/micro structure analysis.
% both kind of files will will look like:
%  P    T    S   sig_t
% ...  ...  ...  ...
%
% usage:  > vmp_p2profiles('Pfiles')
%   or    > vmp_p2profiles('Pfiles', 'gmtd')
%         > vmp_p2profiles('Pfiles', ['fine', 'micro', 'dateformat','clockMismatch',clockOffset])
%            (you can use any combination of terms in brackets )
%
%   varargin signification:
%             - 'fine'          -> create fine*.dat (for reorange use)
%             - 'micro'         -> create micro*.dat (for reorange use)
%             - 'dateformat'    -> outfile in profile_yyyymmdd_001.mat style
%             - 'clockMismatch' -> correct for time offset if needed. The
%                                  time offset (clockOffset), expressed in days, must be 
%                                  entered after the 'clockMismatch' option.
%
%
% ** use view_profile.m to have an idea of the work done...
%

% Author: Frédéric Cyr - 2009/09/03
%
%   - modified by F. Cyr - 2009/12/17:
%       used the modification from split_profile.m to flagged data with NaN
%       when portions of a same profile are stuck together
%   - modified by F. Cyr - 2010-01-11:
%       remove NaN flags for P and p
%   - modified by F. Cyr - 2010-03-23:
%       Now the parameter is not the number of .P files, but a file
%       containing .P files names... more flexible if, for
%       instance, .P files are not regular (ex: DAT001, DAT003,
%       DAT004...)
%   - modified by F. Cyr - 2010-03-23:
%       After modifying p2mat to deal with .P files with missing
%       fields (ex: Galbraith's micro_c), this function need to be
%       modified. For the moment, this modification is quite cheap
%       because I only filled micro_c with NaN values..
%       Noticed also that the code is not performant because
%       there's a load in the "for profile" loop... Could change it
%       only by loading outside the loop and give new name to
%       variable inside the loop...
%
%   - modified by C. Hamel 17/05/2010
%       *turn off p2mat call (line 83) because it had allready be proceed
%       *add a if exist('fluoro') and if exist('trans') because my data
%       dind't containt those elements
%
%   - modified by C. Hamel 17/05/2010
%       put a continue at line 95 if indices==1, that means that the
%       function split_profiles.m had not been able to split the profile
%       because the .p file was not made of good values. split_profile then
%       returned indices=1.
%   - modified by C. Hamel 18/05/2010
%       at line 237: this if/else script ensure that even if the time
%       reference, [round(indices(i,2)/FS)]==0, that the first (and not
%       the 0 elemnt witch doesn't exist) element of the time array will
%       be assigned to the variable 'time'.
%
%   -mofdified by F. Cyr - 2010-05-28
%       comments on line 142-144-145
%
%   - modified by F. Cyr - 2010-11-04
%       add new input "varargin" for the user to have the choice of
%       whether or not save fineXXX.dat and microXXX.dat
%
%   - modified by F. Cyr - 2010-11-04:
%       Branch version where the time is corrected for some days in
%       the july 2010 mission. See just before the saving is done!
%       To fo this, a new varargin input as been added
%   - modified by F. cyr - 2011-12-08:
%       - Remove SBC=conductivity_TPcorr(SBC,SBT,P); since correction
%         is now made in p2mat.m
%       - Add treatment of SP, SA, CT (TEOS-10), see p2mat.m
%       - converstion from micro_c to micro_s is now make by gsw
%       toolkit and not anymore by salinity.m. There is also new
%       variables, ct1, ct2, sa1, sa2, sp1, sp2 (conserv_temp,
%       abs_sal, pract_sal) that will be created, only if
%       profile_lat/lon are in CALIB.DAT
%   - modified by F. cyr - 2012-11-05:
%       Put back split_profile for fluoro, removed by C. Hamel
%   - modified by F. cyr - 2013-02-22:
%       introduced new 'vmp_splitProfiles.m'
%   - modified by F. Cyr - 2013-07-26 ********** RENAMED FUNCTION
%       Now named "vmp_p2profiles.m"
%       Add an option to name output file with date format:
%       profile_yyyymmdd_001.mat
%   - modified by D. Bourgault and F. Cyr - 2013-09-05
%       Make it more casual and more general
%   - modified by F. Cyr - 2015-05-08:
%       Now save *.P original file name ("Pfile_archive") in each profile 
%        (needed for .ODF achiving on the OGSL) 


% usage ex: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Some parameters
wMin = 0.4;  % Minimum acceptable value for fall speed (m/s)
wMax = 0.8;  % Maximum acceptable value for fall speed (m/s)

% Varargin test
if isempty(varargin)==1
    sm            = 0; %save micro = no!
    sf            = 0;
    dateformat    = 0;
    clockMismatch = 0;
elseif size(varargin,2) == 1
    sm = strcmp(varargin{1}, 'micro');
    sf = strcmp(varargin{1}, 'fine');
    clockMismatch = strcmp(varargin{1}, 'clockMismatch');
    dateformat = strcmp(varargin{1}, 'dateformat');
    
elseif size(varargin,2) == 2
    sm = strcmp(varargin{1}, 'micro') + strcmp(varargin{2}, 'micro');
    sf = strcmp(varargin{1}, 'fine') + strcmp(varargin{2}, 'fine');
    clockMismatch = strcmp(varargin{1}, 'clockMismatch') + strcmp(varargin{2}, 'clockMismatch');
    dateformat = strcmp(varargin{1}, 'dateformat') + strcmp(varargin{2}, 'dateformat');
    
elseif size(varargin,2) == 3
    sm = strcmp(varargin{1}, 'micro') + strcmp(varargin{2}, 'micro') ...
       + strcmp(varargin{3}, 'micro');
    sf = strcmp(varargin{1}, 'fine') + strcmp(varargin{2}, 'fine') ...
       + strcmp(varargin{3}, 'fine');
    clockMismatch = strcmp(varargin{1}, 'clockMismatch') + strcmp(varargin{2}, 'clockMismatch') ...
                  + strcmp(varargin{3}, 'clockMismatch');
    dateformat = strcmp(varargin{1}, 'dateformat') + strcmp(varargin{2}, 'dateformat') ...
               + strcmp(varargin{3}, 'dateformat');
elseif size(varargin,2)==4
    sm = strcmp(varargin{1}, 'micro') + strcmp(varargin{2}, 'micro') ...
       + strcmp(varargin{3}, 'micro') + strcmp(varargin{4}, 'micro');
    sf = strcmp(varargin{1}, 'fine') + strcmp(varargin{2}, 'fine') ...
       + strcmp(varargin{3}, 'fine') + strcmp(varargin{4}, 'fine');
    clockMismatch = strcmp(varargin{1}, 'clockMismatch') + strcmp(varargin{2}, 'clockMismatch') ...
                  + strcmp(varargin{3}, 'clockMismatch') + strcmp(varargin{4}, 'clockMismatch');
    dateformat = strcmp(varargin{1}, 'dateformat') + strcmp(varargin{2}, 'dateformat') ...
               + strcmp(varargin{3}, 'dateformat') + strcmp(varargin{4}, 'dateformat');
else
    disp('Wrong input... try "help vmp_p2profiles"')
    return
end

% Find the position of the 'clockMismatch' argument within the varargin cell.
% The following argument is the actual closk offset (in days).
if clockMismatch
    iClockMismatch = find(strcmp(varargin, 'clockMismatch') == 1);
    clockOffset    = varargin{iClockMismatch+1};
end

% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
Pfiles = char(C{1});

no_files = size(Pfiles,1);

profile = 1; %compute the number of the recorded profile (initialisation)

% Save the original working directory in variable ORIG_PATH
ORIG_PATH = pwd;

for count = 1:no_files %loop on DAT000, DAT001, etc.
    fname = Pfiles(count, :);
    I = find(fname == ' ');
    fname(I)         = [];
    fname(end-1:end) = [];
    
    % Get current folder and DAT* folder if different
    slashPos = findstr(fname,'/');
    if isempty(slashPos) == 1
        MY_PATH = './';
        FNAME   = fname;
    else
        MY_PATH = fname(1:slashPos(end));
        FNAME   = fname(slashPos(end)+1:end);
    end
    command = sprintf('cd %s', MY_PATH);
    eval(command);
    
    % Track .P file name (for OGSL archiving)
    Pfile_archive = [FNAME '.P'];
    
    % use p2mat
    vmp_p2mat(FNAME);
    
    % split into profiles
    indices = vmp_splitProfiles(FNAME);
    
    % modified C. Hamel
    if isnan(indices) == 1
        disp('continue')
        continue
    end
    
    no_profiles = length(indices(:,1)); % number of profiles
    
    %eval(['cd ' ORIG_PATH]);
    eval(['cd ' '''',ORIG_PATH,''''])
    
    for i = 1:no_profiles
        
        load(fname)
        
        % Correct time for Fred's mistake during July 2010 mission
        % (see log book)
        if clockMismatch == 1
            MTIME = MTIME + clockOffset;
            mtime = mtime + clockOffset;
            timen = datenum(time(:,1:8)) + clockOffset;
            time(:,1:8) = datestr(timen,13);
            if i == 1
                disp(['Date and time corrected for the following clock mismatch: ',...
                      num2str(clockOffset),' days']);
            end    
        end
        
        if count == 1 % check which date (dateformat option)
            mtime0 = round(mtime(1));
            mdate = datestr(mtime(1), 29);
        end
        
        % Here's the modification to deal with missing fields it is
        % possible to add any other fieldif needed
        if exist('c1_nonphys') == 0 %c1 missing
            c1_nonphys = p.*0;
            c1_nonphys(1:length(c1_nonphys))=NaN;
        end
        if exist('fluoro') == 0 %fluoro missing (C. Hamel)
            fluoro = p.*0;
            fluoro(1:length(fluoro))=NaN;
            fluoro_units = 'none';
        end
        if exist('trans') == 0 %trans missing (C. Hamel)
            trans = p.*0;
            trans(1:length(trans))=NaN;
            trans_units = 'none';
        end
        
        if indices(i,4)==1 %FLAG!
            
            %for the correction we impose NaN where speed doesnt satisfies
            %conditions On some variable,
            
            % We first split W
            W = W(indices(i,2):indices(i,3));
            I = find(W > wMax | W < wMin); % Remove speed which are not betw. 0.4-0.8 m/s
            
            % Split slow variables
            P   = P(indices(i,2):indices(i,3));
            SBC = SBC(indices(i,2):indices(i,3));
            SBS = SBS(indices(i,2):indices(i,3));
            SBT = SBT(indices(i,2):indices(i,3));
            
            % and correction
            SBC(I) = NaN;
            SBS(I) = NaN;
            SBT(I) = NaN;
            MTIME  = MTIME(indices(i,2):indices(i,3));
            
            if exist('profile_lat')~=0
                SP = SP(indices(i,2):indices(i,3)); % TEOS-10
                SA = SA(indices(i,2):indices(i,3)); % TEOS-10
                CT = CT(indices(i,2):indices(i,3)); % TEOS-10
                SP(I) = NaN;
                SA(I) = NaN;
                CT(I) = NaN;
            end
            
            % Split fast variables, same correction
            imin = (indices(i,2)-1)*round(fs/FS)+1;
            imax = indices(i,3)*round(fs/FS);
            
            w = w(imin:imax);
            I = find(w > wMax | w < wMin); % Remove speed which are not betw. 0.4-0.8 m/s
            
            az         = az(imin:imax);
            fluoro     = fluoro(imin:imax);
            p          = p(imin:imax);
            pitch      = pitch(imin:imax);
            roll       = roll(imin:imax);
            shear1     = shear1(imin:imax);
            shear2     = shear2(imin:imax);
            t1_nonphys = t1_nonphys(imin:imax);
            t2_nonphys = t2_nonphys(imin:imax);
            c1_nonphys = c1_nonphys(imin:imax);
            trans      = trans(imin:imax);
            mtime      = mtime(imin:imax);
            
            %az(I)=NaN;
            fluoro(I) = NaN;
            %p(I)=NaN; % Dewey doesnt deal with that
            %pitch(I)=NaN;
            %roll(I)=NaN;
            shear1(I)     = NaN;
            shear2(I)     = NaN;
            t1_nonphys(I) = NaN;
            t2_nonphys(I) = NaN;
            c1_nonphys(I) = NaN;
            trans(I)      = NaN;
            
        else % no corection
            
            % Split slow variables
            P     = P(indices(i,2):indices(i,3));
            SBC   = SBC(indices(i,2):indices(i,3));
            SBS   = SBS(indices(i,2):indices(i,3));
            SBT   = SBT(indices(i,2):indices(i,3));
            W     = W(indices(i,2):indices(i,3));
            MTIME = MTIME(indices(i,2):indices(i,3));
            
            if exist('profile_lat')~=0
                SP = SP(indices(i,2):indices(i,3)); % TEOS-10
                SA = SA(indices(i,2):indices(i,3)); % TEOS-10
                CT = CT(indices(i,2):indices(i,3)); % TEOS-10
            end
            
            % Split fast variables
            imin = (indices(i,2)-1)*round(fs/FS)+1;
            imax = indices(i,3)*round(fs/FS);
            
            az         = az(imin:imax);
            fluoro     = fluoro(imin:imax);
            p          = p(imin:imax);
            pitch      = pitch(imin:imax);
            roll       = roll(imin:imax);
            shear1     = shear1(imin:imax);
            shear2     = shear2(imin:imax);
            t1_nonphys = t1_nonphys(imin:imax);
            t2_nonphys = t2_nonphys(imin:imax);
            c1_nonphys = c1_nonphys(imin:imax);
            trans      = trans(imin:imax);
            w          = w(imin:imax);
            mtime      = mtime(imin:imax);
            
        end
        
        % Calibrate T and C
        t1 = calibrate_micro(t1_nonphys, SBT, W, fs, FS, 0.26, 1); % the distance betw. sensors must be verified
        t2 = calibrate_micro(t2_nonphys, SBT, W, fs, FS, 0.26, 1);
        c  = calibrate_micro(c1_nonphys, SBC, W, fs, FS, 0.38, 1);
        t1_units = SBT_units;
        t2_units = SBT_units;
        
        % computing t as the mean of t1 and t2
        t = mean([t1 t2], 2);
        
        %computing salinity from micro-conductivity
        %s = salinity(p, t, c);
        s = gsw_SP_from_C(c, t, p);
        s_units = 'psu';
        
        % compute sa, sp, ct from GSW toolkit
        if exist('profile_lat') ~= 0 &  exist('profile_lon') ~= 0
            % pratical salinity (SP), Absolute salinity (SA) and conserv. T (CT)
            sp = s;
            [sa, in_ocean] = gsw_SA_from_SP(sp,p,profile_lon,profile_lat);
            ct = gsw_CT_from_t(sa, t, p);
            % units
            sp_units = 'psu'; sa_units = 'g/kg'; ct_units = 'degC';
            clear in_ocean profile_lon profile_lat
        end
        
        % Split super-slow variables
        if (round(indices(i,2)/FS))==0      % modified C.Hamel
            time = time(1:round(indices(i,3)/FS), :);
            date = date(1:round(indices(i,3)/FS), :);
        else
            time = time(round(indices(i,2)/FS):round(indices(i,3)/FS), :);
            date = date(round(indices(i,2)/FS):round(indices(i,3)/FS), :);
        end
        
        
        %save .MAT and .dat (profileXXX.mat, profileXXX.dat)
        if dateformat
            if round(mtime(1)) ~= mtime0;
                profile = 1;
                mtime0  = round(mtime(1));
                mdate   = datestr(mtime(1), 29);
            end
            
            data_fname  = sprintf(['profile_' mdate '_%3.3i.mat'], profile);
            micro_fname = sprintf(['micro_'   mdate '_%3.3i.dat'], profile);
            fine_fname  = sprintf(['fine_'    mdate '_%3.3i.dat'], profile);
            disp(['saving ' data_fname])
            
            
        else % original format profileXXX.mat
            
            data_fname  = sprintf('profile%3.3i', profile);
            micro_fname = sprintf('micro%3.3i.dat', profile);
            fine_fname  = sprintf('fine%3.3i.dat', profile);
            disp(sprintf('Saving profile%3.3i', profile));
            
        end

        if exist('SA') ~= 0
            save(data_fname,'date','time','MTIME','mtime',...  %saving .MAT
                'fs','FS','fs_units','FS_units',...
                'P','p','P_units','p_units',...
                'W','w','W_units','w_units',...
                'pitch','roll','az','pitch_units','roll_units','az_units',...
                'shear1','shear2','shear1_units','shear2_units',...
                't1','t2','t1_units','t2_units', 's', 's_units',...
                'SBT','SBC','SBS','SBT_units','SBC_units','SBS_units',...
                'fluoro','trans','fluoro_units','trans_units', ...
                'SP', 'SA', 'CT', 'SP_units', 'SA_units', 'CT_units', ...
                'sp', 'sa', 'ct', 'sp_units', 'sa_units', 'ct_units', ...
                 'Pfile_archive');
        else
            save(data_fname,'date','time','MTIME','mtime',...  %saving .MAT
                'fs','FS','fs_units','FS_units',...
                'P','p','P_units','p_units',...
                'W','w','W_units','w_units',...
                'pitch','roll','az','pitch_units','roll_units','az_units',...
                'shear1','shear2','shear1_units','shear2_units',...
                't1','t2','t1_units','t2_units', 's', 's_units',...
                'SBT','SBC','SBS','SBT_units','SBC_units','SBS_units',...
                'fluoro','trans','fluoro_units','trans_units', 'Pfile_archive');
        end
        
        % Extcract 4 fields for the reorange method (saved in .dat)
        % Now only if wanted! (nov. 2010)
        if sm==1
            % -- Computing density and sigma-t for microscale-- %
            dens = sw_dens(s, t, p); %from CSIRO
            sig_t = dens-1000;
            micro_fields = [p t s sig_t];
            save(micro_fname, '-ascii', 'micro_fields');  %saving .dat
        end
        if sf==1
            % -- Computing density and sigma-t for finescale-- %
            DENS = sw_dens(SBS, SBT, P);
            SIG_T = DENS-1000;
            fine_fields = [P, SBT, SBS, SIG_T];
            save(fine_fname, '-ascii', 'fine_fields');  %saving .dat
        end
        
        
        clear date time mtime MTIME ...
            fs FS fs_units FS_units ...
            P p P_units p_units ...
            W w W_units w_units ...
            pitch roll az pitch_units roll_units az_units ...
            shear1 shear2 shear1_units shear2_units ...
            t1 t2 t1_units t2_units  s  s_units ...
            SBT SBC SBS SBT_units SBC_units SBS_units ...
            fluoro trans fluoro_units trans_units c1_nonphys
        
        % The following could be adapted to save with a fixed number of digits
        %fid = fopen('profile005.dat', 'w');
        %fprintf(fid, '%3.3f %3.3f %3.3f %3.3f\n', profile005';
        %fclose(fid);
        
        profile = profile + 1; %Increment the name of the next file
        
    end
end
toc
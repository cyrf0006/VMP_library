function vmp_rename_Pfiles(file_names);

%function vmp_rename_Pfiles_and_profiles(file_names);

% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
Pfiles = char(C{1});

no_files = size(Pfiles,1);

profile = 1; %compute the number of the recorded profile (initialisation)

% Save the original working directory in variable ORIG_PATH
ORIG_PATH = pwd;
!rm tmp
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
    
    theDate = fname(slashPos(end-2)+1:slashPos(end-1)-1);
    Pfile_archive = ['DAT_' theDate '_' FNAME(4:end) '.P'];
    CALIB_archive = ['CALIB_' theDate '.DAT'];
    SETUP_archive = ['SETUP_' theDate '.TXT'];
    TEXT_archive = ['DAT_' theDate '.TXT'];

    Pfile_orig = [fname '.P'];
    TEXT_orig = [MY_PATH 'DAT.TXT'];
    SETUP_orig = [MY_PATH 'SETUP.TXT'];
    CALIB_orig = [MY_PATH 'CALIB.DAT'];
    
    % copy .P file
% $$$     command = sprintf('!cp %s %s', Pfile_orig, Pfile_archive);
% $$$     eval(command);
    
    if count == 1
        % copy DAT.TXT file
        fid = fopen(TEXT_orig);
        if fid ~= -1        
            command = sprintf('!cp %s %s', TEXT_orig, TEXT_archive);
            eval(command);
        end
        
        % copy SETUP.TXT and CALIB.DAT
        fid = fopen(SETUP_orig);
        if fid ~= -1
            command = sprintf('!cp %s %s', SETUP_orig, SETUP_archive);
            eval(command);
            command = sprintf('!cp %s %s', CALIB_orig, CALIB_archive);
            eval(command);
        end
    end
    
    command = sprintf('cd %s', MY_PATH);
    eval(command);
    
    % Track .P file name (for OGSL archiving)
    Pfile = [FNAME '.P'];
    
    MATfile = [FNAME '.mat'];
    % load DAT
    fid = fopen(MATfile);
    if fid ~= -1
        load(FNAME);
    else
        disp([fname ' does not exist, CONTINUE'])
        continue
    end
    
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
  
    disp([FNAME sprintf('   --> %d profiles', no_profiles)])
    
    for i = 1:no_profiles        
        command = sprintf('!echo "%s" >> tmp', Pfile_archive);
        eval(command);
    end
end
eval(['cd ' '''',ORIG_PATH,''''])

!mv tmp Pfile_archive.txt


keyboard
        
      
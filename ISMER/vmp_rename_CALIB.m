function vmp_rename_CALIB(file_names);

%function vmp_rename_Pfiles_and_profiles(file_names);

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
    else
        MY_PATH = fname(1:slashPos(end));
    end
    
    theDate = fname(slashPos(end-2)+1:slashPos(end-1)-1);
    CALIB_archive = ['CALIB_' theDate '.DAT'];
    SETUP_archive = ['SETUP_' theDate '.TXT'];
    TEXT_archive = ['DAT_' theDate '.TXT'];

    TEXT_orig = [MY_PATH 'DAT.TXT'];
    SETUP_orig = [MY_PATH 'SETUP.TXT'];
    CALIB_orig = [MY_PATH 'CALIB.DAT'];
    
    
    command = sprintf('!cp %s %s', TEXT_orig, TEXT_archive);  
    eval(command);
    command = sprintf('!cp %s %s', SETUP_orig, SETUP_archive);
    eval(command);
    command = sprintf('!cp %s %s', CALIB_orig, CALIB_archive);
    eval(command);

    
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
    
    eval(['cd ' '''',ORIG_PATH,''''])
  
    disp([SETUP_archive])

end
eval(['cd ' '''',ORIG_PATH,''''])



keyboard
        
      
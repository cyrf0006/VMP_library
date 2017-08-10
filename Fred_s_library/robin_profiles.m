% to be executed in ~/WINDEX/data_processing

disp('Profiles...')
file_names = 'profile_list/rikiProfiles.list';
epsfile_names = 'profile_list/rikiEpsProfiles.list';


fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
Pfiles = char(C{1});
noFiles = size(Pfiles,1);

fid = fopen(epsfile_names);
C = textscan(fid, '%s', 'delimiter', '\n');
EPSfiles = char(C{1});


!mkdir robin_profiles
for i = 1:noFiles

    disp(sprintf('no. %d / %d', i, noFiles))
    
    if i<10
        outFile = sprintf('profile000%d.mat', i);
        outFile2 = sprintf('epsprofile000%d.mat', i);
    elseif i<100
        outFile = sprintf('profile00%d.mat', i);
        outFile2 = sprintf('epsprofile00%d.mat', i);
    elseif i<1000
        outFile = sprintf('profile0%d.mat', i);
        outFile2 = sprintf('epsprofile0%d.mat', i);
    else
        outFile = sprintf('profile%d.mat', i);
        outFile2 = sprintf('epsprofile%d.mat', i);
    end
   
    % profiles
    fname = Pfiles(i, :);
    I = find(fname==' ');   
    fname(I) = [];  
    command1 = sprintf('cp %s ./robin_profiles/%s', fname, outFile);
    
    % epsprofiles
    fname = EPSfiles(i, :);
    I = find(fname==' ');   
    fname(I) = [];  
    command2 = sprintf('cp %s ./robin_profiles/%s', fname, outFile2);
    
    
    %disp(command1);
    system(command1);
    system(command2);
end

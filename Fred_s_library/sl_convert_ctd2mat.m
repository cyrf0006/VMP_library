function sl_convert_ctd2mat(ctdFiles)

% Loop for converting .cnv Files to .mat, using cnv2mat.m 
% usage ex: convert_ctd2mat('ctdfiles_20130609')
% where 'ctdfiles_20130609' is from (ex:):
% 'ls -1 SBE19plus_01906786_2013_06_10*.cnv > ctdfiles_20130610'

% F. Cyr - june 2013
    
fid = fopen(ctdFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

noProfiles = size(files,1);

for i = 1:noProfiles
    
    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    OUTFILE = [fname(1:end-4) '.mat'];

    disp(sprintf('%s -> %s', fname, OUTFILE));

    [lat,lon,mtime,data,names,sensors]=cnv2mat(fname);

    save(OUTFILE, 'lat', 'lon', 'mtime', 'data', 'names', 'sensors')
end

    
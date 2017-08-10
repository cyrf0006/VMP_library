function profile_count(file_names)


% load *.P files names (file in which are recorded .P files)
fid = fopen(file_names);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

no_profile = size(files,1);

totalDepth = 0;

%%%%%%%%%%%%%%%%%%%%
% loop on profiles %
%%%%%%%%%%%%%%%%%%%%

for i = 1:no_profile

    if mod(i, 100) == 0;
        disp(sprintf('%d out of %d', i, no_profile));
    end
    
    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    load(fname)
    totalDepth = totalDepth+max(p_eps1);
end



disp(sprintf('Total depth: %d km in %d profiles', totalDepth/1000, no_profile))
    %%%%%%%%%%%%%%%%%%%%%%
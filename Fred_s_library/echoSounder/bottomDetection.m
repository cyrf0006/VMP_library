function bottomDetection(infiles, mindepth, outfile)
    

fid = fopen(infiles);
C = textscan(fid, '%s', 'delimiter', '\n');
dt4files = char(C{1});

noFiles = size(dt4files, 1);


bottom = [];
mtime = [];
for i = 1:noFiles
    
    fname = dt4files(i, :);
    I = find(fname==' ');   
    fname(I) = [];    
    
    dtx=rddtx(fname);
    
    J = find(dtx.range>mindepth);
    
    
    
    [Y, I] = max(dtx.vals(J,:));
    bot = [];
    for i = 1:length(I)          
        bot(i) = dtx.range(J(I(i))); 
    end

    bottom = [bottom bot];
    mtime = [mtime dtx.imtime];
    
    
end

% Attempt to remove unreal bottom
L = 100;
for i = 1:length(bottom)
    I = i-L/2:i+L/2;
    J = find(I<=0 | I>length(bottom));
    I(J) = [];
    M = nanmean(bottom(I));
    if bottom(i)>1.05*M | bottom(i)<.95*M
        bottom(i) = M;
    end
    
end



bottomTime = mtime;
save(outfile, 'bottom', 'bottomTime'); 

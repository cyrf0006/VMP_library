z = 2:293;


for i = 1:length(z)
    
    [Y, I] = min(abs(z(i)-P));
    timesec(i) = (MTIME(I)-MTIME(z(1)))*86400;
end

mm = floor(timesec./60);
ss = timesec - 60*mm;

[z' round(timesec)' mm' round(ss)']

    
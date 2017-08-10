function time2 = time2hightide(timeVec, tidefile)


% usage ex:
%    time2 = time2hightide(timeVec, '/home/cyrf0006/WINDEX/data_processing/BMix_study/tide_2009-2012.dat');
%
    
tide  = load(tidefile);

% time 2 hightide
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);

% find high tide time 
clear T L
count = 1;
for i = 2:length(mtime)-1

    if level(i)>level(i-1) & level(i)>=level(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        count = count+1;
    end
end
% T contains the hour of each high tide
clear mtime

time2 = nan(size(timeVec));
for i = 1:length(timeVec)
    [Y, I] = min(abs(timeVec(i)-T));
    time2(i) = (timeVec(i)-T(I))*24;
end

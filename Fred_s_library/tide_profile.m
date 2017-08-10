function tide_profile(tidefile, prof_files)
    
% Returns when profiles have been taken relative to the closest
% high tide. T is thus a vector (mtime style) containing how far
% from the closest highest tide.
% ex: tide_profile('tide_2009-2011.dat', 'file_names_veryall')   
% ex2: tide_profile('tide_2009-2011.dat', 'file_names_riki') 
%                            or
%      tide_profile('tide_2009-2011.dat', 'file_names_bndry') 

springtide = 3.5; %m
    

tide  = load(tidefile);
    
mtime = datenum(tide(:,1), tide(:,2), tide(:,3), tide(:,4), tide(:,5), 0);
level = tide(:,6);


% find high tide time 
count = 1;
for i = 2:length(mtime)-1

    if level(i)>level(i-1) & level(i)>level(i+1)
        T(count) = mtime(i); % high tide time
        L(count) = level(i); % high tide level
        count = count+1;
    end
end
% T contains the hour of each high tide
clear mtime

% load *.P files names (file in which are recorded .P files)
fid = fopen(prof_files);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});

count2=1;
count3=1;
for i = 1:size(files,1)
    fname=files(i,:);
    I = find(fname==' ');   
    fname(I) = [];
    load(fname);

    if exist('mtime_eps')==1
        [Y, I] = min(abs(mtime_eps(1)-T));
        A(i) = (mtime_eps(1)-T(I))*24;
        B(i) = L(I); %level of the closest hightide
    else
        [Y, I] = min(abs(mtime(1)-T));
        A(i) = (mtime(1)-T(I))*24;
        B(i) = L(I); %level of the closest hightide
    end
    
    if L(I) > springtide
 
        spring_files(count2,:) = files(i,:);
        count2 = count2+1;
    else
        neap_files(count3,:) = files(i,:);
        count3 = count3+1;
    end
    
end

figure(1)
hist(A, 3)
save 'time_to_hightide.mat' A

% neap-spring
figure(2)
hist(B, 3)
title('max tide during sampling period')

figure(3)
hist(L, 3)
title('max tide of the full neap-spring cycle')

%keyboard

dlmwrite('springfiles', spring_files, 'delimiter', '');
dlmwrite('neapfiles', neap_files, 'delimiter', '');



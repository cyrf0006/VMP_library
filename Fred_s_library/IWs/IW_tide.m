function IW_tide(tidefile, IW_file)
  
% function IW_tide(tidefile, IW_file)
%
% ex: IW_tide('../tide_shear/tide_2009-2011.dat', 'IW_detected.mat')
%
%

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


load(IW_file)

keyboard

count2=1;
count3=1;
for i = 1:size(IW,2)
    [Y, I] = min(abs(T-IW(2,i)));
    A(i) = (IW(2,i)-T(I))*24;
    B(i) = L(I); %level of the closest hightide
    if L(I) > springtide 
        spring(count2) = i;
        count2 = count2+1;
    else
        neap(count3) = i;
        count3 = count3+1;
    end
end



figure(1)
hist(A, 13)
title('Sept. 12 - Oct. 19')
xlabel('time to hightide')
ylabel('"IW" occurence')
print('-dpng', '-r300','IW_histo_all.png')


% neap-spring
figure(2)
hist(A(spring), 13)
title('spring tide')
xlabel('time to hightide')
ylabel('"IW" occurence')
print('-dpng', '-r300','IW_histo_spring.png')


figure(3)
hist(A(neap), 13)
title('neap tide')
xlabel('time to hightide')
ylabel('"IW" occurence')
print('-dpng', '-r300','IW_histo_neap.png')



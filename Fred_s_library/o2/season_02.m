
clear
month = 4:11;

hyp = 2;



for i = 1:length(month) 
    
    
    % -- Load and store bootstraped climatology -- %
    if month(i) <10
            o2filename  = sprintf('O2_bootclim_0%d.dat', month(i));
     else
            o2filename  = sprintf('O2_bootclim_%d.dat', month(i));
     end
    
    o2profile = load(o2filename);

    
    if i == 1
        depth = o2profile(:,1);
        dz = depth(2)-depth(1);
    end
    
    n(i) = datenum(999, month(i), 15, 0,0,0); % climatology, 15th of the month, year 999
    O2mat(i,:) = o2profile(:,2);
    O2_97p2mat(i,:) = o2profile(:,4); % 95% CI
    O2_2p5mat(i,:) = o2profile(:,3);

    % find Hypoxic layer
    I = find(o2profile(:,2)<=hyp);
    hypox(i) = nanmean(o2profile(I,2));
    hypox_2p5(i) = nanmean(o2profile(I,3));
    hypox_97p5(i) = nanmean(o2profile(I,4));
    
    % find layer >200m
    I = find(depth>=200);
    o200(i) = nanmean(o2profile(I,2));
    o200_2p5(i) = nanmean(o2profile(I,3));
    o200_97p5(i) = nanmean(o2profile(I,4));
    
    
end



figure(1)
plot(n, hypox, 'k', 'linewidth', 2)
hold on 
plot(n, hypox_2p5, '--k', 'linewidth', 0.5)
plot(n, hypox_97p5, '--k', 'linewidth', 0.5)
datetick('x', 3) 
ylabel('[O_2] (mg/l)')
title('mean concentration when hypoxic')


figure(2)
plot(n, o200, 'k', 'linewidth', 2)
hold on 
plot(n, o200_2p5, '--k', 'linewidth', 0.5)
plot(n, o200_97p5, '--k', 'linewidth', 0.5)
datetick('x', 3) 
ylabel('[O_2] (mg/l)')
title('mean concentration under 200m')
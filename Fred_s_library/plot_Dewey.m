basic = 0;

if basic
    figure(2)
    clf
    set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 15 12])
    
    loglog(k1, ps1, 'k', 'linewidth', 2)
    hold on
    plot(uk, uphi, '--g', 'linewidth', 2);  
    
    plot([kl0 kl0], [1e-6 1e-1], '--r');
    plot([ku0 ku0], [1e-6 1e-1], '--r');  
    xlim([1e-1 1e3])
    ylim([1e-6 1e-1])
    text(2e-1, 5e-2, sprintf('eps0 = %0.2g', eps0_fc))
    text(2e-1, 3e-2, sprintf('eps1 = %0.2g', eps1_fc))
    text(2e-1, 2e-2, sprintf('eps2 = %0.2g', eps2_fc))
    %text(2e-1, 1.2e-2, sprintf('epsNa = %d', 7.5*viscosity*uvar))
    
% $$$ plot([uk(uil) uk(uil)], [1e-6 1e-1], '--m');  
% $$$ plot([uk(uiu) uk(uiu)], [1e-6 1e-1], '--m');  

    title(sprintf('z = %1.1f; L1 = %d; L2 = %d; L3 = %d', zfft(i), iloop, ie, in))


    disp('press any key to continue')
    pause
else % subplots
    if iloop == 1 & ie == 1 & in == 1
        figure(2)
        clf
        set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 20 22])
   
        count = 0;
        % *********************** Adjust_space.m ************************ %
        % Fields required by the function adjust_space.m. Please fill every
        % of the following and call "adjust_space" in the script whenever
        % you want. Do not touch four last fields
        ncol = 2; % no. subplot column
        nrow = 3; % no. subplot row
        dx = 0.02 ; % horiz. space between subplots
        dy = 0.03; % vert. space between subplots
        lefs = 0.08; % very left of figure
        rigs = 0.01; % very right of figure
        tops = 0.04; % top of figure
        bots = 0.05; % bottom of figure
        figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
        figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
        count_col = 1;
        count_row = 1;
        % *************************************************************** %
    end
    
    isplot=1;
    if iloop == 1 & ie == 1 & in == 2
        s1 = subplot(3,2,1);
    elseif iloop == 1 & ie == 2 & in == 2
        s2 = subplot(3,2,2);
    elseif iloop == 1 & ie == 3 & in == 2
        s3 = subplot(3,2,3);
    elseif iloop == 1 & ie == 4 & in == 2
        s4 = subplot(3,2,4);
    elseif iloop == 2 & ie == 1 & in == 2
        s5 = subplot(3,2,5);
    elseif iloop == 3 & ie == 1 & in == 2
        s6 = subplot(3,2,6);
    else
        isplot=0;
    end

    if isplot
           
        loglog(k1, ps1, 'k', 'linewidth', 2)
        hold on
        plot(uk, uphi, '--g', 'linewidth', 2);  
        
        plot([kl0 kl0], [1e-10 1e-1], '--r');
        plot([ku0 ku0], [1e-10 1e-1], '--r');  

        % for 6.5m
        if zfft(i) == 6.5
            text(2e-1, 5e-2, sprintf('eps0 = %0.2g', eps0_fc))
            text(2e-1, 3e-2, sprintf('eps1 = %0.2g', eps1_fc))
            text(2e-1, 2e-2, sprintf('eps2 = %0.2g', eps2_fc))
            title(sprintf('L1 = %d; L2 = %d; L3 = %d', iloop, ie, in))

            xlim([1.5e-1 5e2])
            ylim([5e-6 1e-1])
        elseif zfft(i) == 54.5
            % for 54.5m
            text(2e-1, 5e-6, sprintf('eps0 = %0.2g', eps0_fc))
            text(2e-1, 3e-6, sprintf('eps1 = %0.2g', eps1_fc))
            text(2e-1, 2e-6, sprintf('eps2 = %0.2g', eps2_fc))
            title(sprintf('L1 = %d; L2 = %d; L3 = %d', iloop, ie, in))

            xlim([1.5e-1 5e2])        
            ylim([3e-9 5e-5])            
        end
        
        adjust_space
        count = count+1;

        if count == 1
            ylabel('\phi (s^{-2} cpm^{-1})');
        elseif count == 5
            xlabel('k(cpm)');
        end

        if count == 6;
            
            set(s1,'xticklabel', []);
            set(s2,'xticklabel', []);
            set(s3,'xticklabel', []);
            set(s4,'xticklabel', []);    
            set(s2,'yticklabel', []);
            set(s4,'yticklabel', []);
            set(s6,'yticklabel', []);
            
            print('-dpng', 'Dewey_algo.png')
            set(gcf, 'renderer', 'painters')
            print('-depsc2', 'Dewey_algo.eps')
        end

    end
end

    
    
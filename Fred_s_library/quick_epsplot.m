function quick_epsplot(no_profile)

% for visual inspection of every epsilon profile

for profile=1:no_profile
    
    if profile<10
        eps_fname = sprintf('eps_profile00%d', profile);
    else
        if profile<100
            eps_fname = sprintf('eps_profile0%d', profile);
        else %profile>100
            eps_fname = sprintf('eps_profile%d', profile);
        end
    end
    
    
    load(eps_fname)
    plot(eps1, p_eps1)
    set(gca, 'ydir', 'reverse')
    hold on
    plot(eps2, p_eps2, 'r')
    hold off
    
    pause
    
end

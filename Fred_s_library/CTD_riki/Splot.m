clear

leg = ['b', 'r', 'k', 'g', 'm', 'y', 'c'];
month = [5 6 7 8 9 10 11];

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 10])

for i = 1:length(month)
    
    if month(i)<10
        Tfname  = sprintf('T_bootclim_0%d.dat', month(i));
        Sfname  = sprintf('S_bootclim_0%d.dat', month(i));
    else
        Tfname  = sprintf('T_bootclim_%d.dat', month(i));
        Sfname  = sprintf('S_bootclim_%d.dat', month(i));
    end
    
    S = load(Sfname);
    plot(S(:,2), S(:,1), leg(i));
    if i==1
       set(gca, 'ydir', 'reverse')
       xlabel('S')
       ylabel('P(dbars)')
       hold on
    end
     
end


figure(1)
hold off
legend('May', 'jun', 'Jul', 'aug', 'Sept', 'oct', 'Nov', 'location', 'southwest')

print(gcf, '-depsc2', 'Splot_may-nov.eps')
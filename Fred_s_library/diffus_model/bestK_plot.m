clf


bestK1 = load('bestK_T_mean.dat');
bestK2 = load('bestK_T_plusmin.dat');
bestK3 = load('bestK_T_minplus.dat');

[Y,I] = min(bestK2(:,2));
m1 = bestK2(I,1);
[Y,I] = min(bestK3(:,2));
m2 = bestK3(I,1);
[Y,I] = min(bestK1(:,2));
mm = bestK3(I,1);

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 10 10])

plot(bestK1(:,1), bestK1(:,2), 'k', 'linewidth', 2)
hold on
plot(bestK2(:,1), bestK2(:,2), 'k', 'linewidth', 2)
plot(bestK3(:,1), bestK3(:,2), 'k', 'linewidth', 2)

% error region
plot([m1 m1], [0 350], '--k', 'linewidth', 1)
plot([m2 m2], [0 350], '--k', 'linewidth', 1)
plot([mm mm], [0 250], '--k', 'linewidth', 1)

hold off

xlabel('K (m^2 s^{-1})')
ylabel('misfit')
axis([2e-5 9e-5 0 400])



print('-deps2', 'bestK_error.eps')



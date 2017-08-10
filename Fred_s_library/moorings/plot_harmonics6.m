% Tf
yy = [1e-5 1e-2];
xx = 1./[15.981/24  15.981/24];
loglog(xx, yy, 'k')
text(xx(1), 2e-2, 'f', 'color', 'k', 'horizontalalignment', 'center')



%M2
yy = [1e-8 1e-1];
xx = 1./[12.42/24  12.42/24];
loglog(xx, yy, 'k')
text(xx(1), 1e-9, 'M_2', 'color', 'k', 'horizontalalignment', 'center')

%M4
yy = [1e-5 1e-1];
xx = 1./([12.42/24  12.42/24]./2);
loglog(xx, yy, 'k')
text(xx(1), 2e-1, 'M_4', 'color', 'k', 'horizontalalignment', 'center')

% K1
yy = [1e-8 1e-1];
xx = 1./([23.93/24 23.93/24]);
loglog(xx, yy, 'k')
text(xx(1), 1e-9, 'K_1', 'color', 'k', 'horizontalalignment', 'center')


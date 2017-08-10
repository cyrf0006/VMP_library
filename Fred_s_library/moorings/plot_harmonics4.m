% Tf
yy = [1e-5 1e-2];
xx = 1./[15.981/24  15.981/24];
loglog(xx, yy, 'k')
text(xx(1), 2e-2, 'T_f', 'color', 'k', 'horizontalalignment', 'center')

% Tf/2
yy = [1e-5 1e-2];
%yy = [1e0 1e3];
xx = 1./([15.981/24  15.981/24]/2);
loglog(xx, yy, 'k')
text(xx(1), 2e-2, 'T_f/2', 'color', 'k', 'horizontalalignment', 'center')



%M2
yy = [1e-5 1e-1];
%yy = [1e0 1e4];
xx = 1./[12.42/24  12.42/24];
loglog(xx, yy, 'k')
text(xx(1), 2e-1, 'M_2', 'color', 'k', 'horizontalalignment', 'center')

%M3
yy = [1e-5 1e-2];
%yy = [1e0 1e4];
%xx = 1./([12.42/24  12.42/24]*2/3);
xx = 1./[12.42/24  12.42/24] + 1./[23.93/24  23.93/24];
loglog(xx, yy, 'k')
text(xx(1), 5e-6, 'M_2 + K_1', 'color', 'k', 'horizontalalignment', 'center')


%M4
yy = [1e-5 1e-1];
%yy = [1e0 5e3];
xx = 1./([12.42/24  12.42/24]./2);
loglog(xx, yy, 'k')
text(xx(1), 2e-1, 'M_4', 'color', 'k', 'horizontalalignment', 'center')


%M5
yy = [1e-5 5e-2];
%yy = [1e0 5e3];
xx = 1./[12.42/24  12.42/24] + 2.*(1./[23.93/24  23.93/24]);
loglog(xx, yy, 'k')
text(xx(1), 5e-6, 'M_2 + 2K_1', 'color', 'k', 'horizontalalignment', 'center')


%M6
yy = [1e-5 1e-1];
%yy = [5e-1 5e2];
xx = 1./([12.42/24  12.42/24]./3);
loglog(xx, yy, 'k')
text(xx(1), 2e-1, 'M_6', 'color', 'k', 'horizontalalignment', 'center')

% K1
yy = [1e-5 1e-1];
%yy = [5e-1 5e2];
xx = 1./([23.93/24 23.93/24]);
loglog(xx, yy, 'k')
text(xx(1), 5e-6, 'K_1', 'color', 'k', 'horizontalalignment', 'center')

% $$$ % Max N
% $$$ yy = [1e-7 1e-3];
% $$$ %yy = [5e-1 5e2];
% $$$ xx = [800 800]
% $$$ loglog(xx, yy, '--k')
% $$$ %text(xx(1), 2e-1, 'mean ', 'color', 'k', 'horizontalalignment', 'center')
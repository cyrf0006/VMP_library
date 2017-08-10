% $$$ % --------------------------------------------- %
% $$$ % -------------- DIFFUSE_MODEL.m -------------- %
% $$$ % --------------------------------------------- %
% $$$ 
% $$$ %diffus_model('T_climato_04.dat', 'S_climato_04.dat' , datenum(999,4,22), 9e-7, 'constant_e', 'daily_heatflx_clim.dat', 100, 3600*24*30, 3600*24*30*6)
% $$$ 
% $$$ diffus_model('T_climato_04.dat', 'S_climato_04.dat' ,datenum(999,4,22), 5e-5, 'constant_k', 'daily_heatflx_clim.dat',100, 3600*24, 3600*24*30*7)
% $$$ 
% $$$ load approx_e.dat
% $$$ K=approx_e;
% $$$ 
% $$$ %diffus_model('T_climato_04.dat', 'S_climato_04.dat' , datenum(999,4,22), [K(:,1) K(:,2)], 'variable_e', 'daily_heatflx_clim.dat', 100, 3600*24*30, 3600*24*30*6)
% $$$ 
% $$$ 
% $$$ load meanK_IML4_5m-bin
% $$$ 
% $$$ %diffus_model('T_climato_04.dat', 'S_climato_04.dat' ,
% $$$ %datenum(999,4,22), [P_K_bin K_mean_bin], 'variable_k',
% $$$ %'daily_heatflx_clim.dat', 100, 3600*24, 3600*24*30*7)



% ---------------------------------------------------- %
% -------------- DIFFUSE_MODEL_noflux.m -------------- %
% ---------------------------------------------------- %

%diffus_model_noflux('T_bootclim_04.dat', 'S_bootclim_04.dat' ,datenum(999,4,15), 5.5e-5, 'constant_k', 'daily_forcing.dat',100, 3600*24, 3600*24*30*8)

% $$$ 
% $$$ %KK = load('K_all_boot.dat');
% $$$ %KK = load('K_veryall.dat');
% $$$ %KK = load('K_veryveryall.dat');
% $$$ KK = load('K_veryall_corr.dat');
% $$$ K = [KK(:,1), KK(:,2)];
% $$$ %K = [KK(:,1), KK(:,3)];
% $$$ diffus_model_noflux('T_bootclim_04.dat', 'S_bootclim_04.dat', datenum(999,4,15), [K(:,1) K(:,2)], 'variable_k','daily_forcing.dat', 100, 3600*24, 3600*24*30*8)


% ---------------------------------------------------- %
% -------------- DIFFUSE_MODEL_O2.m -------------- %
% ---------------------------------------------------- %

% $$$ %diffus_model('T_climato_04.dat', 'S_climato_04.dat' , datenum(999,4,22), 9e-7, 'constant_e', 'daily_heatflx_clim.dat', 100, 3600*24*30, 3600*24*30*6)
% $$$ 

load SOD 
t1 = datenum(999, 1, 1, 0, 0, 0);
tt = t1+t./86400;
dlmwrite('SOD.dat', [tt' SOD'],'delimiter',' ','precision',12)


% $$$ % ---- test a ----- %
% $$$ P = [200:300]';
% $$$ O2 = ones(101,1).*120;
% $$$ dlmwrite('O2_init.dat', [P O2],'delimiter',' ','precision',12)
% $$$ 
% $$$ 
% $$$ diffus_model_O2('O2_init.dat' ,datenum(999,1,1, 0, 0, 0), 1e-4, 'constant_k', 'SOD.dat',900, 3600*24, 3600*24*810, 120)
% $$$ 


% ---- test b ---- %
% $$$ P = [150:300]';
% $$$ O2 = ones(length(P),1).*120;
% $$$ I = find(P>=150 & P<=200);
% $$$ O2(I) = -2*P(I)+520;
% $$$ 
% $$$ dlmwrite('O2_init.dat', [P flipud(O2)],'delimiter',' ','precision',12)
% $$$ diffus_model_O2('O2_init.dat' ,datenum(999,1,1, 0, 0, 0), 1e-4, 'constant_k', 'SOD.dat',900, 3600*24, 3600*24*810, 220)

% $$$ % ---- test c ---- %
% $$$ P = [150:300]';
% $$$ O2 = ones(length(P),1).*120;
% $$$ I = find(P>=150 & P<=200);
% $$$ O2(I) = -2*P(I)+520;
% $$$ 
% $$$ dlmwrite('O2_init.dat', [P flipud(O2)],'delimiter',' ','precision',12)
% $$$ diffus_model_O2('O2_init.dat' ,datenum(999,1,1, 0, 0, 0), 1e-5, 'constant_k', 'SOD.dat',900, 3600*24, 3600*24*810, 220)



% ---- test Lefort ---- %
P = [150:525]';
O2 = ones(length(P),1).*170;
%I = find(P>=150 & P<=200);
%O2(I) = -2*P(I)+520;
dlmwrite('O2_init_lefort.dat', [P flipud(O2)],'delimiter',' ','precision',12)

load SOD.dat
A = SOD(:,2);
A = A*0+3540;
A = A*1000/86400/365;
SOD(:,2) = A;
dlmwrite('SOD_lefort.dat', SOD,'delimiter',' ','precision',12)


diffus_model_O2_lefort('O2_init_lefort.dat' ,datenum(999,1,1, 0, 0, 0), 1e-4, 'constant_k', 'SOD_lefort.dat',900, 3600*24, 3600*24*810, 170)




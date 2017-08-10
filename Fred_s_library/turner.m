% in /home/cyrf0006/PhD/CTD_IML4/BIN/renamed

load S_bootclim_07.dat
load T_bootclim_07.dat
S = S_bootclim_07(:,2)
T = T_bootclim_07(:,2)
P = S_bootclim_07(:,1);


SA = gsw_SA_from_SP(S, P, -70, 48);
CT = gsw_CT_from_t(SA,T, P);

% $$$ alpha = gsw_alpha(SA,CT,P);
% $$$ beta = gsw_beta(SA,CT,P);
% $$$ g = 9.81;
% $$$ 
% $$$ dTdz = gradient(CT, 1);
% $$$ dSdz = gradient(SA, 1);
% $$$ 
% $$$ 
% $$$ N_theta = g.*alpha.*dTdz;
% $$$ N_S = g.*beta.*dSdz;
% $$$ Tu = 

[Tu, Rsubrho, p_mid] = gsw_Turner_Rsubrho(SA,CT,P);
plot(Tu, P(1:end-1)+.5)
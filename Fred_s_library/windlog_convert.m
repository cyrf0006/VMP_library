function U = windlog_convert(U0, Z1, Z2);

% This function convert wind from a determined height to another one using
% a log wind profile
% usage ex: wind10 = windlog_convert(wind2, 2, 10);
% to convert wind at 2m to wind at 10m
%
% Author: F. Cyr - march 2010
%
% ---------------------------------------------------------------------

if Z1 ==0 | Z2 ==0
    disp('Z1 or Z2 can''t be zero!')
    return
end


% parameters
kappa = 0.41;
cd_black = [1.46;1.7;1.8;1.6;1.5]; % for U10 = [15-20; 20-22; 22-25; 25-30; >30]
cd_pond = [1.2; 0.49; 0.065]; %[1.2; 0.49 + 0.065U10] for U10 = [0-11; 11-15];

U = U0;
C_d = U0.*0; % initialization

%1st iteration
    
% Drag coef calculation;
I = find(U<11);
C_d(I) = cd_pond(1);
I = find(U>=11 & U < 15);
C_d(I) = cd_pond(2) + cd_pond(3).*U(I);
I = find(U>=15 & U < 20);
C_d(I) = cd_black(1);
I = find(U>=20 & U < 22);
C_d(I) = cd_black(2);
I = find(U>=22 & U < 25);
C_d(I) = cd_black(3);
I = find(U>=25 & U < 30);
C_d(I) = cd_black(4);
I = find(U>=30);
C_d(I) = cd_black(5);

C_d = C_d./1000; % in 1e-3


F = 1-sqrt(C_d)./kappa.*log(Z2/Z1); % factor betw. 2 wind speed (U(Z1)/U(Z2))
U = U0./F; %new wind speed at Z2
U_prec = U0;

%other iterations
while std(U-U_prec)>0.1
    
    U_prec=U;
    
    % Drag coef calculation;
    I = find(U<11);
    C_d(I) = cd_pond(1);
    I = find(U>=11 & U < 15);
    C_d(I) = cd_pond(2) + cd_pond(3).*U(I);
    I = find(U>=15 & U < 20);
    C_d(I) = cd_black(1);
    I = find(U>=20 & U < 22);
    C_d(I) = cd_black(2);
    I = find(U>=22 & U < 25);
    C_d(I) = cd_black(3);
    I = find(U>=25 & U < 30);
    C_d(I) = cd_black(4);
    I = find(U>=30);
    C_d(I) = cd_black(5);
    
    C_d = C_d./1000; % in 1e-3

    F = 1-sqrt(C_d)./kappa.*log(Z2/Z1); % factor betw. 2 wind speed (U(Z1)/U(Z2))

    U = U0./F; %new wind speed at Z2

end

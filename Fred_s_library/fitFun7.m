function out = fitFun7(a)    
% optimizing function
    
global Z E

[Y, I] = sort(Z);
out = sum((E(I) - f(a, Y)).^2);

%out = sum((E - f(a, Z)).^2);

 
function e = f(a, z)
    
    e = a(1)./z + a(2);
    


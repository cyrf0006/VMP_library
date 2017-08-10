function out = fitFun4(a)    
% optimizing function
    
global Z U Ub
   
out = sum((U - f(a, Ub, Z)).^2);

 
function u = f(a, ub, z)
    
    kappa = 0.41;    
    u = a(1)./kappa.*log(z./a(2)) + a(2)*a(1).^2/1e-6;


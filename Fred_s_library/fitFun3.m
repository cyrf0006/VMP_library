function out = fitFun3(a)    
% optimizing function
    
global Z U Ub
   
out = sum((U - f(a, Ub, Z)).^2);

 
function u = f(a, ub, z)
    
    kappa = 0.41;    
    u = sqrt(a(1))./kappa.*ub.*log(z./a(2));
    


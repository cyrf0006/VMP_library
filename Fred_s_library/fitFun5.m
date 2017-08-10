function out = fitFun5(a)    
% optimizing function
    
global Z E
   
out = sum((E - f(a, Z)).^2);

 
function e = f(a, z)
    
    kappa = 0.41;    
    e = a.^3/kappa./z;
    


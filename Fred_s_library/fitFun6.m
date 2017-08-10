function out = fitFun6(a)    
% optimizing function
    
global Z E
   
out = sum((E - f(a, Z)).^2);

 
function e = f(a, z)
    
    e = a(1)*z + a(2);
    


function out = fitFun9(a)    
% optimizing function
    
global Z T
   
out = sum((T - f(a, Z)).^2);

 
function t = f(a, z)

    t = a(1)*z.^2 + a(2)*z + a(3);
    


function out = fitFun8(a)    
% optimizing function
    
global Z T
   
out = sum((T - f(a, Z)).^2);

 
function t = f(a, z)
    
    t = a(1)*z.^5 + a(2)*z.^4 + a(3)*z.^3 + a(4)*z.^2 + a(5)*z + a(6);
    


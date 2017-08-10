function out = fitFun10(a)    
% optimizing function
    
global Z S

out = sum((S - f(a, Z)).^2);

 
function s = f(a, z)
    
    s = a(1) + a(2).*tanh((z-a(3))./a(4));
    


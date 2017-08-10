function out = fitFun2(a)    
% optimizing function
    
    global Z S
   
    out = sum((S - f(a, Z)).^2);
end

 
function s = f(a, z)
    s = a(1).*exp(a(2)./(z+a(3)));
end
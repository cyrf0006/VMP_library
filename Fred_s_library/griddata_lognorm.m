function [X, Y, Z] = griddata_lognorm(x, y, z, xVec, yVec);

% usage ex: 
%     [XI,YI,Z] = griddata_lognorm(LN2,LS2,LEPS,X,Y);

    
dx = xVec(2) - xVec(1);
dy = yVec(2) - yVec(1);
 
[X,Y] = meshgrid(xVec, yVec);
Z = nan(size(X));

for i = 1:length(xVec)    
    for j = 1:length(yVec)
        I = find(log10(x)>=xVec(i)-dx/2 & log10(x)<xVec(i)+dx/2 & log10(y)>=yVec(j)-dy/2 & log10(y)<yVec(j)+dy/2); 

        if ~isempty(I) == 1
            m = nanmean(log(z(I)));
            s2 = sum((log(z(I))-m).^2)/length(I);
            Z(j,i) = exp(m+s2/2);
            %Z(j,i) = nanmean(z(I));            
        end        
    end
end

        
Z = log10(Z);
        
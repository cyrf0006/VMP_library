function build_M(matfile, cond)
%%% --------------------------------------------- %%%
% This function will create a matrix containing the dissipation rate for
% each overturn found

%clear
load(matfile)

depth = 2:0.5:180;
depth = single(depth); % 2 digits precision for the following find
tide = -6:0.1:6;
tide = single(tide); % 2 digits precision for the following find

dim1 = length(depth); % 1801
dim2 = length(tide); % 121

if cond == 0 %Condition to see if there's a existing file
    mat = sparse(dim1, dim2); %Here mat will be a matrix with dissipation over depth and tide cycle
else
    load M
end

no_ovt = length(dissip_matrix(:,1));


for ovt = 1:no_ovt
    
    
    zmin = single(round(mean(dissip_matrix(ovt, 3)*10))); %(round every 10cm)
    zmax = single(round(mean(dissip_matrix(ovt, 4)*10))); %(round every 10cm)    
    
    t = single(round(mean(dissip_matrix(ovt, 2)*10))); %(round every 0.1hr)
    
    %OVT = sparse(dim1, 1);
%    depthi1 = find(depth==zmin/10); %indice of the beginning of the overturn
%    depthi2 = find(depth==zmax/10); %indice of the end of the overturn   
    
    depthi1 = dsearchn(depth', zmin/10); %indice of the beginning of the overturn
    depthi2 = dsearchn(depth', zmax/10); %indice of the beginning of the overturn   
    
    %OVT(depthi1:depthi2) = dissip_matrix(ovt, 1);
    OVT(1:depthi2-depthi1+1) = dissip_matrix(ovt, 1);
    
    timei = find(tide==t/10);
  
    mat(depthi1:depthi2,timei) = OVT;
    
    clear OVT zmin zmax t depthi1 depthi2
    
end

save M mat 
function dissip_ovt(no_profile, tide_file)
% dissip_ovt(no_profile, tide_file)
% 
% This function creates a "dissipation matrix" containing information from
% many profiles. This dissipation matrix will look like:
% [ epsilon  time-vs-HT  Zmin  Zmax ]
%     ...        ...      ...   ...  
%
% where, 
% - epsilon: TKE dissipation
% - time-vs-HT: the moment relative to high tide
% - Zmin: beginning of the overturn
% - Zmax: end of the overturn
%
% In input, the function needs the number of profile to analyse and a
% tide_file containing the hour, minutes and level of the highest high tide
% of the day:
%
% usage ex: 
% dissip_ovt(9, 'HT_09-16')
% where 9 is the number of profile and HT_09-16 contains a vector called 
% hht = [HH MM LEVEL] (hour / minutes / level in m)
%
% The script create a .m file with a name relative to tide_file, in this
% case it would be: dissip_matrix_09-16.mat
%
% Author: Frederic Cyr - 2009/01/12
%
%
%% -------------------------------------------------------------- %%

% % some constants
nu = 1e-6;
g = 9.81;%m^2/s
rho_0 = 1.035e3;%kg/m^3

load(tide_file);
count_ovt = 1;

for profile = 1:no_profile

    % the name of the profile
    if profile<10
        data_fname = sprintf('profile00%d', profile);
        eps_fname = sprintf('eps_profile00%d', profile);
        fine_ovt = sprintf('ovt.fine00%d.ovt', profile);        
        disp(sprintf('profile00%d', profile));
    else
        if profile<100
            data_fname = sprintf('profile0%d', profile);
            eps_fname = sprintf('eps_profile0%d', profile);
            fine_ovt = sprintf('ovt.fine0%d.ovt', profile);   
            disp(sprintf('profile0%d', profile));
        else %profile>100
            data_fname = sprintf('profile%d', profile);
            eps_fname = sprintf('eps_profile%d', profile);
            fine_ovt = sprintf('ovt.fine%d.ovt', profile);            
            disp(sprintf('profile%d', profile));
        end
    end
    
    
    % We load the consider profile and overturn
    load(data_fname);
    load(fine_ovt);
        
    %compute the number of ovt for the consider profile
    ovt_size = size(ovt);
    no_ovt = ovt_size(1); 
    
    %loop on overturns
    for o = 1:no_ovt
        zmin = ovt(o, 3);
        zmax = ovt(o, 4);
                
        Imid = nearestpoint((zmax+zmin)/2, P); %indice of the closest/middle point of the ovt
        I1 = nearestpoint(zmin, P); % coWINDEX/data_processing/2009_07_21uld be change by dsearchn.m
        I2 = nearestpoint(zmax, P);
        time_ind = round(Imid/64); %because time is recorded every 1s
                
        if time_ind~=0 %make sure that the overturn is not at surface !!!!!! DOUBLE CHECK THIS!
             
        time_ovt = time(time_ind, :); %the mean time of the ovt (string vector)
        t = [str2num(time_ovt(1, 1:2)) str2num(time_ovt(1, 4:5)) str2num(time_ovt(1, 7:8))]; % [hh mm ss] (num vector)
        
        % -- Dissipation -- %
        load(eps_fname);
        
        p_eps = mean([p_eps1; p_eps2]); %pas ideal si pas meme longueur ou decalee (devrait faire un test si un est 0 on prend l<autre)
        eps = mean([eps1; eps2]);
        ind = dsearchn(p_eps', (zmax+zmin)/2); %indice of the dissipation matrix p_eps corresponding to the needed depth
        epsilon = eps(ind); %dissipation corresponding to this indice
        
                       
        % -- Time vs high tide -- %
        if t(3)>29
            t(2)=t(2)+1; %round to minutes
        end
        
        time_ovt_dec = t(1)+t(2)/60; %time in decimal
        time_hht_dec = hht(1)+hht(2)/60; 
                
        ovt_vs_hht = time_ovt_dec - time_hht_dec;
        
        % transform time in a 12-hour cycle ([-6, 6])
        if ovt_vs_hht > 6
            ovt_vs_hht = ovt_vs_hht - 12;
        elseif ovt_vs_hht < -6
            ovt_vs_hht = ovt_vs_hht + 12;
        end
            
        
        % -- making dissip_matrix -- %
        dissip_matrix(count_ovt, :) = [epsilon ovt_vs_hht zmin zmax];
        count_ovt = count_ovt+1;
        
        % -- Computing Ozmidoz and Thorpe scales for comparison -- %
        
        L_t(count_ovt) =  ovt(o, 5);  % thorpe scale already computed in reorange
        
        ind1 = ovt(o, 1); % indices of the overturn
        ind2 = ovt(o, 2);       
            
        DENS = sw_dens(SBS, SBT, P); % density profile
        N2 = -g/rho_0*(DENS(ind2)-DENS(ind1))/(zmax-zmin); % Buoy. freq. as a linear fit wihtin the overturn
        L_o(count_ovt) = sqrt(epsilon/(N2)^(3/2));
        
        end
        
    
    end
end

      out_file = ['dissip_matrix_' tide_file(4:8)];
      save (out_file, 'dissip_matrix') 
      

        
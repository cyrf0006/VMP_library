    % monthly climato interpolated to daily values
    clear
    time_model = load('N_diffus_daily.dat');
    time_clim = datenum(999, 1:12, 15);
    
    Tclim = nan(300,12);
    for month = 4:11
         if month <10
            Toutname  = sprintf('T_climato_0%d.dat', month);
         else
            Toutname  = sprintf('T_climato_%d.dat', month);
         end

         Tmonth = load(Toutname);
         Tclim(:,month)=Tmonth(:,2);
    end
    
    for i = 1:size(Tclim,1)
        
        Tdaily(i,:) = interp1(time_clim, Tclim(i,:), time_model);
    end
    
    
    dlmwrite('T_diffus_daily_tobs.dat', Tdaily,'delimiter',' ','precision',6)      

function view_eps_profiles(no)

% function view_profile(no)
%
% For a quick view at allo TKE dissipation profiles generated!
%
% Where 
% - no: the number of profile to look at
% 
% ex: view_profile(50)
% will plot 50 profiles of eps_1 & eps_2 against p1 & p2 and will ask you to hit any key
% between 2 plot. This function could be use to have a quick look of
% profiles and may be use to detect profile splitting errors...
%
% Author: Frederic Cyr - 2010/05/25
%
% ------------------------------------------------------------------- %


for profile= 1:no
        
       if profile<10
        data_fname = sprintf('eps_profile00%d', profile);
       else
           if profile<100
            data_fname = sprintf('eps_profile0%d', profile);
         else %profile>100
            data_fname = sprintf('eps_profile%d', profile);
         end
       end
    
       load(data_fname);
            
       figure(1)
       clf
       
       semilogx(eps1, p_eps1)
       set(gca, 'ydir', 'reverse');
       hold on
       semilogx(eps2, p_eps2, 'r')
       hold off
       title(sprintf('profile %d',profile))
       disp('Press any key to continue \n');
       pause 

end
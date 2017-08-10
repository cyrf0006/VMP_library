function boot_ctd_cmos(ctd_files)

% "ls -1 *.bin > ctd_files"




% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. 
ncol = 4; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.05 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])
    
% -- Some parameters -- %
dz = 1; % would need to be adjusted if not already binned
P = [1:dz:300]';   
nboot = 500;
rho_0 = 1.025e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
Tcil = 1; %degC, threshold for the CIL core
freeze_pt = -1.8;

% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

no_files = size(ctd,1); 

% compute time vector for considered profiles
% Only work for Peter's Perl renaming script
n = datenum(str2num(ctd(:,1:4)), str2num(ctd(:,6:7)), str2num(ctd(:,9:10)),str2num(ctd(:,12:13)), str2num(ctd(:,15:16)), 0 );

% Monthly average for which months?
%month = unique(str2num(datestr(n,5)));
month = 4:11;


for i = 1:length(month)
    disp(sprintf('month %d', month(i)))
    II = find(str2num(datestr(n,5))==month(i));

    clear T S
    for j = 1:length(II)       
        
        file = load(ctd(II(j),:));
        
        % Quality control + save into matrix T and S
        if ~isempty(file)==1 % check if file empty
            
            pres = file(:,1);
            temp = file(:,2);
            sali = file(:,3);
            
            %remove -99
            I = find(temp==-99);
            temp(I) = NaN;
            I = find(sali==-99);
            sali(I) = NaN;
                
            if max(pres) >= max(P) & min(pres)==min(P)
                I = find(file(:,1)>=min(P) & file(:,1)<=max(P));
                
                % fill matrix with empty value
                T(P,j)=NaN;
                S(P,j)=NaN;
                
                T(pres(I),j) = temp(I);
                S(pres(I),j) = sali(I);               

            else % CTD profile doesnt start at p=1m or pmax smaller than 300m
                
                I = find(file(:,1)>=max(min(pres), min(P)) & file(:,1)<=min(max(pres), max(P)));
                % fill matrix with empty value
                T(P,j)=NaN;
                S(P,j)=NaN;
                
                % substituate empty value where it is possible
                T(pres(I),j) = temp(I);
                S(pres(I),j) = sali(I);
            
            end
                 
        else %file empty
            T(P,j)=NaN;
            S(P,j)=NaN;
        end
        
    end  % for j
    
    % -- Quality control -- %
    % Remove empty columns
    I = find(sum(~isnan(T),1)>0);
    T = T(:,I);
    I = find(sum(~isnan(S),1)>0);
    S = S(:,I);
    
    % vertical interpolation to remove NaNs
    for j = 1:size(T,2)
        I = find(isnan(T(:,j))==1);
        if ~isempty(I)==1 
            p = P;
            t = T(:,j);
            t(I) = [];
            p(I) = [];
            T(:,j) = interp1(p, t, P);
        end
    end
    
    for j = 1:size(S,2)
        I = find(isnan(S(:,j))==1);
        if ~isempty(I)==1
            p = P;
            s = S(:,j);
            s(I) = [];
            p(I) = [];
            S(:,j) = interp1(p, s, P);
        end
    end
    
    % -- bootstrap -- %
    disp('bootstrap...')
    NT = size(T,2); % NT and NS may not be equal...
    NS = size(S,2);
    for b = 1:nboot
        
        % temp
        r = rand(NT,1);
        r = ceil(r*NT/1);
        T_boot_b(:,b) = nanmean(T(:,r),2);

        % salinity
        r = rand(NS,1);
        r = ceil(r*NS/1);
        S_boot_b(:,b) = nanmean(S(:,r),2);
    end
   
    % --- mean values --- %
    %T, S
    T_boot_dot = nanmean(T_boot_b,2);
    S_boot_dot = nanmean(S_boot_b,2);    
    
    % Erosion variable (mean)
    [Y, J] = min(T_boot_dot);
    cilT_boot_dot = Y;
    cilZ_boot_dot = P(J);
    CIL = find(T_boot_dot<=Tcil);
    if  ~isempty(CIL) == 1
        DENS = sw_dens(S_boot_dot(CIL), T_boot_dot(CIL), P(CIL)); 
        cilH_boot_dot =  cp*nansum(DENS.*(T_boot_dot(CIL)-freeze_pt))*dz./length(CIL);
        cilD_boot_dot = length(CIL);  
    else
        cilH_boot_dot = NaN;
        cilD_boot_dot = 0;
    end
        
    % --- Compute error --- %
    %T,S
    T_error = sqrt(sum((diff([T_boot_b T_boot_dot],1, 2)).^2, 2)./(nboot-1));
    S_error = sqrt(sum((diff([S_boot_b S_boot_dot],1, 2)).^2, 2)./(nboot-1));

    % Erosion variable (error)
    [Y, J] = min(T_boot_b);
    cilT_boot_b = Y;
    cilZ_boot_b = J;
    for k = 1:nboot
        CIL = find(T_boot_b(:,k)<=Tcil);
        if  ~isempty(CIL) == 1
            DENS = sw_dens(S_boot_b(CIL, k), T_boot_b(CIL, k), P(CIL)); 
            cilH_boot_b(k) =  cp*nansum(DENS.*(T_boot_b(CIL, k)-freeze_pt))*dz./length(CIL);
            cilD_boot_b(k) = length(CIL);
        else
            cilH_boot_b(k) = NaN;
            cilD_boot_b(k) = 0;
        end   
    end
  
    cilT_error(month(i)) = sqrt(sum((diff([cilT_boot_b cilT_boot_dot],1, 2)).^2, 2)./(nboot-1));
    cilZ_error(month(i)) = sqrt(sum((diff([cilZ_boot_b cilZ_boot_dot],1, 2)).^2, 2)./(nboot-1));
    cilH_error(month(i)) = sqrt(sum((diff([cilH_boot_b cilH_boot_dot],1, 2)).^2, 2)./(nboot-1));
    cilD_error(month(i)) = sqrt(sum((diff([cilD_boot_b cilD_boot_dot],1, 2)).^2, 2)./(nboot-1));

    % --- Compute 95% confidence interval --- %
    T_sort = sort(T_boot_b, 2);
    S_sort  = sort(S_boot_b, 2); 
   
    CI_2p5 = round(2.5/100*nboot);
    CI_97p5 = round(97.5/100*nboot);

    %T,S
    T_2p5 = T_sort(:,CI_2p5);
    S_2p5 = S_sort(:,CI_2p5);
    T_97p5 = T_sort(:,CI_97p5);
    S_97p5 = S_sort(:,CI_97p5);
    

    % Erosion variables (2.5%, 97.5%)   
     for k = 1:nboot         
         [Y, J] = min(T_boot_b(:,k));
         cilT_boot_b(k) = Y;
         cilZ_boot_b(k) = P(J);
         CIL = find(T_boot_b(:,k)<=Tcil);
        if  ~isempty(CIL) == 1
            DENS = sw_dens(S_boot_b(CIL, k), T_boot_b(CIL, k), P(CIL)); 
            cilH_boot_b(k) =  cp*nansum(DENS.*(T_boot_b(CIL, k)-freeze_pt))*dz./length(CIL);
            cilD_boot_b(k) = length(CIL);
        else
            cilH_boot_b(k) = NaN;
            cilD_boot_b(k) = 0;
        end   
    end
  
    cilT_sort = sort(cilT_boot_b, 2);
    cilZ_sort = sort(cilZ_boot_b, 2);
    cilH_sort = sort(cilH_boot_b, 2);
    cilD_sort = sort(cilD_boot_b, 2);

    cilT_2p5(month(i)) = cilT_sort(CI_2p5);
    cilZ_2p5(month(i)) = cilZ_sort(CI_2p5);
    cilH_2p5(month(i)) = cilH_sort(CI_2p5);
    cilD_2p5(month(i)) = cilD_sort(CI_2p5);
    cilT_97p5(month(i)) = cilT_sort(CI_97p5);
    cilZ_97p5(month(i)) = cilZ_sort(CI_97p5);
    cilH_97p5(month(i)) = cilH_sort(CI_97p5);
    cilD_97p5(month(i)) = cilD_sort(CI_97p5);

    T_ave = T_boot_dot;
    S_ave = S_boot_dot;
    cilT_ave(month(i)) = cilT_boot_dot;
    cilZ_ave(month(i)) = cilZ_boot_dot;  
    cilH_ave(month(i)) = cilH_boot_dot; 
    cilD_ave(month(i)) = cilD_boot_dot; 
    
    
    % temperature plot
    figure(1)
    subplot(2, 4, i)
    plot(T_ave, P, 'k', 'linewidth', 1)
    hold on
    x1 = T_97p5;
    x2 = T_2p5;
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(T_ave, P, 'k', 'linewidth', 1)

    %identify cIL
    I = find(T_ave<1);
    if ~isempty(I)==1
        x1 = T_ave(I);
        x2 = T_ave(I)*0+1;
        patch([x1; flipud(x2); x1(1)], [P(I); flipud(P(I)); P(I(1))], ...
              [    0.5212    0.5212    0.7905], 'edgecolor', 'none');
    end
    
    hold off
    set(gca, 'ydir', 'reverse')
    set(gca, 'ygrid', 'on')
    set(gca, 'xgrid', 'on')
    set(gca, 'box', 'on')
    axis([-1 6 0 300])
    ylabel('')
    xlabel('')
    %title(datestr(n(II(j)), 3))
    text(-0.5, 295,datestr(n(II(j)), 3), 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', ...
                        'left','BackgroundColor',[1 1 1])

    if i == 1 | i == 5
        ylabel('P (dbar)')
    end
    if i == 2 | i == 3 | i == 4 | i == 6 | i == 7 | i == 8
        set(gca, 'yticklabel', '')
    end
    if i == 1 | i == 2 | i == 3 | i == 4
        set(gca, 'xticklabel', '')
    end
     if i == 5
        xlabel('T(^{\circ}C)')
    end
    %adjust_space  
    
    
    % salinity plot
    figure(2)
    subplot(2, 4, i)
    plot(S_ave, P, 'k', 'linewidth', 1)
    hold on    
    x1 = S_97p5;
    x2 = S_2p5;
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(S_ave, P, 'k', 'linewidth', 0.5)
    hold off
    set(gca, 'ydir', 'reverse')
    set(gca, 'ygrid', 'on')
    set(gca, 'xgrid', 'on')
    set(gca, 'box', 'on')
    axis([25 35 0 300])
    ylabel('')
    xlabel('')
    %title(datestr(n(II(j)), 3))
    text(26, 295,datestr(n(II(j)), 3), 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', ...
                        'left','BackgroundColor',[1 1 1])
    
    if i == 1 | i == 5
        ylabel('P (dbar)')
    end
    if i == 2 | i == 3 | i == 4 | i == 6 | i == 7 | i == 8
        set(gca, 'yticklabel', '')
    end
    if i == 1 | i == 2 | i == 3 | i == 4
        set(gca, 'xticklabel', '')
    end
    if i == 5
        xlabel('S')
    end
    
    % Save monthly profile
    if month(i) <10
        Toutname  = sprintf('T_bootclim_0%d.dat', month(i));
        Soutname  = sprintf('S_bootclim_0%d.dat', month(i));
    else
        Toutname  = sprintf('T_bootclim_%d.dat', month(i));
        Soutname  = sprintf('S_bootclim_%d.dat', month(i));
    end

    adjust_space    
    
    
    
    
end %for i


figure(1)
%print('-deps2', 'monthly_climT_cmos.eps')
print('-dpng', '-r300','monthly_climT_cmos.png')


figure(2)
print('-dpng', '-r300', 'monthly_climS_cmos.png')


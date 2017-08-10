function boot_o2(ctd_files)

% function boot_ctd(ctd_files)
%
% This function performs bootstrap analysis and monthly average of
% profiles in a list. It also computes the confidence intervals on
% the monthly average. There is 3 main features:
% 
% 1) it saves monthly average with the confidence intervals in
% files like: T_bootclim_0%d.dat and % S_bootclim_0%d.dat
% 2) Plot monthly profiles with 95%CI and save them in
% monthly_climT.eps and monthly_climS.eps.
% 3) Finally it also saves CIL_erosionCI.dat and CIL_erosionCI.dat
% where erosions variables are store with their 95CI. 
%
% input: ctd_files, the list of CTD profiles to use to build
% monthly averages. To obtain them:
%
% NOTE: adjust_sapce.m could not be run on both figure at the same
% time, typically, I run this script twice, and I adjust space only
% one at a time
%
% "ls -1 *.bin > ctd_files"
%
% Author: Frederic Cyr, April 2011
% 
% ------------------------------------------------------------- %

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. 
ncol = 4; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
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

% -- Some parameters -- %
zbin = 1; % would need to be adjusted if not already binned
zmin=1;
zmax=330;
P = [zmin:zbin:zmax]';   
nboot = 500;

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

    clear O2
    for j = 1:length(II)       
        
        file = load(ctd(II(j),:));
        
        % Quality control + save into matrix T and S
        if ~isempty(file)==1 % check if file empty
            
            if size(file, 2)~=5
                continue
            end
            
            pres = file(:,1);
            oxyg = file(:,5);
            
            %remove -99
            I = find(oxyg==-99);
            oxyg(I) = NaN;

            % bin profiles
            for k = 1:length(P)
                I = find(pres>(zmin+(k-1)*zbin) & pres<(zmin+(k)*zbin));
                oxyg_bin(k, :) = nanmean(oxyg(I));
            end
            oxyg = oxyg_bin;
            pres = [round(min(pres)):zbin:round(max(pres))]';
 
            if max(pres) >= max(P) & min(pres)<=min(P)
                I = find(pres>=min(P) & pres<=max(P));
                
                % fill matrix with empty value
                O2(P,j)=NaN;
                
                O2(pres(I),j) = oxyg(pres(I));

            else % CTD profile doesnt start at p=1m or pmax smaller than 300m
                
                %I = find(file(:,1)>=max(min(pres), min(P)) & file(:,1)<=min(max(pres), max(P)));
                I = find(pres>=max(min(pres), min(P)) & pres<=min(max(pres), max(P)));
                % fill matrix with empty value
                O2(P,j)=NaN;
                
                % substituate empty value where it is possible
                O2(pres(I),j) = oxyg(pres(I));
            end
                 
        else %file empty
            O2(P,j)=NaN;
        end
        
    end  % for j
    
    % -- Quality control -- %
    % Remove empty columns
    I = find(sum(~isnan(O2),1)>0);
    O2 = O2(:,I);
    
    % ml/l to mg/l
    O2 = O2.*1.42903;
    
    % Eliminate P > 318m
    I = find(P>318);
    O2(I,:) = nan(zmax-318, size(O2,2));
    
    % vertical interpolation to remove NaNs
    for j = 1:size(O2,2)
        I = find(isnan(O2(:,j))==1);
        if ~isempty(I)==1 
            p = P;
            o = O2(:,j);
            o(I) = [];
            p(I) = [];
            O2(:,j) = interp1(p, o, P);
        end
    end
    
    % -- bootstrap -- %
    disp('bootstrap...')
    NO = size(O2,2); % NT and NS may not be equal...
    for b = 1:nboot        
        r = rand(NO,1);
        r = ceil(r*NO/1);
        O2_boot_b(:,b) = nanmean(O2(:,r),2);
    end
   
    % --- mean values --- %
    O2_boot_dot = nanmean(O2_boot_b,2);
   

    % --- Compute 95% confidence interval --- %
    O2_sort = sort(O2_boot_b, 2);
   
    CI_2p5 = round(2.5/100*nboot);
    CI_97p5 = round(97.5/100*nboot);

    O2_2p5 = O2_sort(:,CI_2p5);
    O2_97p5 = O2_sort(:,CI_97p5);

    O2_ave = O2_boot_dot;
    
    
    
   
    % Oxygen plot
    figure(1)
    subplot(2, 4, i)
    plot(O2_ave, P, 'k', 'linewidth', 1)
    hold on

    x1 = O2_97p5;
    x2 = O2_2p5;
    
    nann = find(isnan(x1)==1);
    if ~isempty(nann)==1 % remove NaNs
    
        x1(nann)=[];
        x2(nann)=[];
        PP = P;
        PP(nann)=[];
        patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.8 .8 .8], 'edgecolor', 'none');
        plot(O2_ave, P, 'k', 'linewidth', 1)
        hold off
    else
        patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
        plot(O2_ave, P, 'k', 'linewidth', 1)
        hold off
    end
    
   
    set(gca, 'ydir', 'reverse')
    set(gca, 'ygrid', 'on')
    set(gca, 'xgrid', 'on')
    set(gca, 'yminortick', 'on')
    set(gca, 'xminortick', 'on')
    set(gca, 'box', 'on')
    axis([0 12 0 320])
    ylabel('')
    xlabel('')
    %title(datestr(n(II(j)), 3))
    text(8, 310,datestr(n(II(j)), 3), 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment',   'left','BackgroundColor',[1 1 1])

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
        xlabel('diss. O_2(mg/l)')
    end
    
    % adjust space!
    adjust_space 
    
    % Save monthly profile
    if month(i) <10
        O2outname  = sprintf('O2_bootclim_0%d.dat', month(i));
    else
        O2outname  = sprintf('O2_bootclim_%d.dat', month(i));
    end

    dlmwrite(O2outname, [P O2_ave O2_2p5 O2_97p5],'delimiter',' ','precision',6)

end %for i



figure(1)
print('-deps2', 'monthly_climO.eps')
% $$$ 
% $$$ 
% $$$ % Save erosion variables (4 variables of fig.6 + error)
% $$$ dlmwrite('CIL_erosion.dat', [cilT_ave; cilT_error; cilD_ave; ...
% $$$                         cilD_error; cilH_ave/1000; cilH_error/1000; cilZ_ave; cilZ_error],'delimiter',' ','precision',6)
% $$$ 
% $$$ % Save erosion variables (4 variables of fig.6 + 2.5% + 97.5%)
% $$$ dlmwrite('CIL_erosionCI.dat', [cilT_ave; cilT_2p5; cilT_97p5; cilD_ave; ...
% $$$                         cilD_2p5; cilD_97p5; cilH_ave/1000; ...
% $$$                     cilH_2p5/1000; cilH_97p5/1000; cilZ_ave; cilZ_2p5; ...
% $$$                    cilZ_97p5],'delimiter',' ','precision',6)
% $$$ 

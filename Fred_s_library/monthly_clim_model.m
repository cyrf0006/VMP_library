function monthly_clim_model(month, model_T, model_N, varargin)

% function monthly_clim_model(month, model_output, varargin)
% 
% Used to plot a monthly climatology of CTD profiles with the
% uncertainty, together with the results of the diffusivity model
% for each monthly average.
%
% usage ex: 
% >> monthly_clim_model([5 6 7 8 9 10 11], 'T_diffus_daily_Kobs.dat', 'N_diffus_daily.dat')
% >> monthly_clim_model([4 5 6 7 8 9 10 11], 'T_diffus_daily_5p5e-5.dat', 'N_diffus_daily.dat', 'T_diffus_daily_Kobs.dat')
% >> monthly_clim_model([4 5 6 7 8 9 10 11],'T_diffus_daily_Kveryveryall.dat', 'N_diffus_daily.dat', 'T_diffus_daily_Kobs.dat')
%
% If only one modeloutput is specified, the results will be a
% subplot for each month where the CTD observations are in a shade
% of gray while the modeled profile if superimposed. If another
% model output is specified, the second profile is also
% superimposed with a dashed line. Note that times of both output
% must be the same... (This is easy to adjust but I dont want to
% include too many parameters for the moment!)
%
% NOTE: to be run in  ~/WINDEX/data_processing/1D_mixing/
%
% Frederic Cyr - march 2011
%
% -------------------------------------------------------------- %

% --- Deal with varargin --- %
if isempty(varargin)==1
    xtra = 0;
elseif size(varargin,2)==1
    xtra = 1;
    model_T2 = varargin{1};
 elseif size(varargin,2)==2
    xtra = 2;
    model_T2 = varargin{1};
    model_T3 = varargin{2};
else
    disp('Wrong input... try "help monthly_clim_model"')
    return
end

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


% -- Zoom options -- %
pmin = 50;
pmax = 150;
pmin = 1;
pmax = 300;

% --- Plot --- %
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

% load model results
Tmod = load(model_T);
Nmod = load(model_N);

% loop on climatology
for i = 1:length(month)
    disp(sprintf('month %d', month(i)))
    II = find(str2num(datestr(Nmod,5))==month(i));
        
    if month(i) <10
        Tname  = sprintf('T_bootclim_0%d.dat', month(i));
    else
        Tname  = sprintf('T_bootclim_%d.dat', month(i));
    end
    
    T = load(Tname);
    Tm = nanmean(Tmod(:,II), 2);
    To_ave = T(:,2);
    To_2p5 = T(:,3);
    To_97p5 = T(:,4);
    P = T(:,1);   

    x1 = To_2p5;
    x2 = To_97p5 ;
    
    % temperature plot
    figure(1)
    subplot(2, 4, i)
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
    hold on
    plot(Tm, P, 'k', 'linewidth', 0.5)
    hold off
    set(gca, 'ydir', 'reverse')
    set(gca, 'ygrid', 'on')
    set(gca, 'xgrid', 'on')
    set(gca, 'box', 'on')
    set(gca, 'xminortick', 'on')
    set(gca, 'yminortick', 'on')

    axis([-1 6 pmin pmax])
    ylabel('')
    xlabel('')
    %    title(datestr(Nmod(II(1)), 3))
    text(-0.5, 295,datestr(Nmod(II(1)), 3), 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', ...
                        'left','BackgroundColor',[1 1 1])

    if i == 1 | i == 5
        ylabel('P (dbar)')
    end
    if i == 5
        xlabel('T(^{\circ}C)')
    end
    if i == 2 | i == 3 | i == 4 | i == 6 | i == 7 | i == 8
        set(gca, 'yticklabel', '')
    end
    if i == 1 | i == 2 | i == 3 | i == 4
        set(gca, 'xticklabel', '')
    end
    
    
    % -- Another model output? -- %
    if xtra == 1 | xtra == 2;
        Tmod2 = load(model_T2);      
        Tm2 = nanmean(Tmod2(:,II), 2);
        hold on
        plot(Tm2, P, '--k', 'linewidth', 0.5)
        hold off
    end
    
    if xtra == 2;
        Tmod3 = load(model_T3);      
        Tm3 = nanmean(Tmod3(:,II), 2);
        hold on
        plot(Tm3, P, '--k', 'linewidth', 0.5)
        hold off
    end
   
    adjust_space
    
% $$$     % salinity plot
% $$$     figure(2)
% $$$     subplot(2, 4, i)
% $$$     plot(S_ave, P, 'k', 'linewidth', 1)
% $$$     hold on
% $$$     plot(S_ave-S_error, P, '--k')
% $$$     plot(S_ave+S_error, P, '--k')
% $$$     hold off
% $$$     set(gca, 'ydir', 'reverse')
% $$$     set(gca, 'ygrid', 'on')
% $$$     set(gca, 'xgrid', 'on')
% $$$     axis([25 35 0 300])
% $$$     ylabel('')
% $$$     xlabel('')
% $$$     title(datestr(n(II(j)), 3))
% $$$ 
% $$$     if i == 1 | i == 5
% $$$         ylabel('P (dbar)')
% $$$     end
% $$$     if i == 1
% $$$         xlabel('S')
% $$$     end
% $$$     if i == 2 | i == 3 | i == 4 | i == 6 | i == 7 | i == 8
% $$$         set(gca, 'yticklabel', '')
% $$$     end
% $$$     if i == 5 | i == 6 | i == 7 | i == 8
% $$$         set(gca, 'xticklabel', '')
% $$$     end
% $$$     
    % Save monthly profile
    
end %for i

figure(1)
print('-deps2', 'monthly_clim_obs_mod.eps')

% $$$ figure(2)
% $$$ print('-deps2', 'monthly_climS.eps')

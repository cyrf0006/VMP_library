function interanual_CIL_stats_2010(ctd_files, year)
    
% New version of interannual_CIL_stats.m. Now a function that take
% in account a list of profiles an years for which plot the
% interanual quantities.
%
% usage ex: interanual_CIL_stats_2010('ctd_files', [1993:2010]')

% Parameters to edit
%ctd_files = 'ctd_files';
%year = [1993:2010]';
month=[4:11]';

mean_type = 4; %  2 = 2weekly; 4 = monthly

TIK = [datenum(999, month,1, 0,0,0)]; %ticklabel for plot
TIK15 = [datenum(999, month,15, 0,0,0)]; %ticklabel for plot
T_LIM = [datenum(999,min(month),15, 0,0,0) datenum(999,max(month),15, 0,0,0)]; %XLIM
                                                              %for
                                                              %plot(apr15
                                                              %-15nov)
T_LIM2 = [datenum(999,min(month),1, 0,0,0) datenum(999,max(month),30, 0,0,0)]; %XLIM for plot
no_days = T_LIM(2)-T_LIM(1);
no_days2 = T_LIM2(2)-T_LIM2(1);
plotlegend = ['b'; 'r'; 'k'; 'y'; 'm'; 'g'; 'c'; 'b']; %for T-S plot


figw = 14; %cm width
figh = 4; %cm height


% Parameters for figure costumization
offset1 = 0.01; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 0.0005; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom
offset5 = 0; % Xtra offset top of figure
offset6 = 0; % Xtra offset right of figure
% $$$ cbar_width = 0.03;
% $$$ clabel_width = 0.08;
% $$$ ti_cbar_frac = 15/16; %reduction of distance bet. colorbar and its
% $$$                       %title
Ylab_offset = -14; %offset between yaxis and ylabel
Ylab_tightoffset = 0.09; %TightInset that must be added for title


% some constant
g = 9.81;
rho_0 = 1.035e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
T_core = 1; %degC, threshold for the CIL core
dz = 1;%m
theta=35;%degrees
freeze_pt = -1.8;

Pres = 1:300; %pressure vector


% ----------------------------------------------------------- % 
% ------------------- Compute time serie -------------------- %
% ----------------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%
% -- CIL properties -- %
%%%%%%%%%%%%%%%%%%%%%%%%

% load all casts

% 2011 version
[T, S, N] = ctd_matrix(ctd_files);

% old version
% $$$ T = load('tprofiles.dat');
% $$$ S = load('sprofiles.dat');
% $$$ N = load('datprofiles.dat');

yyyy = str2num(datestr(N,10));
mm = str2num(datestr(N,5));
dd = str2num(datestr(N,7));


count = 1; % counter for the number of considered pts
    for j=month(1):month(length(month))
for i = year(1):year(length(year))
        

        if mean_type == 2 % 2-weeks average
            
            % first 2 weeks
            I = find(yyyy == i & mm == j & dd<=15);
                
            if ~isempty(I)==1 %nothing to do if I empty
                T_inter(:,count) = nanmean(T(:,I),2);
                S_inter(:,count) = nanmean(S(:,I),2);
                N_inter(count) = datenum(i, j, 7); % 7th day of the month 
            
                count = count+1;
            end
            
            % last 2 weeks
            I = find(yyyy == i & mm == j & dd>15);
            
            if ~isempty(I)==1
                T_inter(:,count) = nanmean(T(:,I),2);
                S_inter(:,count) = nanmean(S(:,I),2);
                N_inter(count) = datenum(i, j, 21); % 21st day of the month 
            
                count = count+1;
            end
            
            
        else % monthly average
            
            I = find(yyyy == i & mm == j);
            
            if ~isempty(I)==1
                T_inter(:,count) = nanmean(T(:,I),2);
                S_inter(:,count) = nanmean(S(:,I),2);
                N_inter(count) = datenum(i, j, 15); % 15th day of the month 
            
                count = count+1;
            end
            
        end
        
    end
end


% ----------------------------------------------------------- % 
% -------------------  Multi-year plots  -------------------- %
% ----------------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%
% -- Tmin warming -- %
%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ figure(2)
% $$$ clf
% $$$ figure(3)
% $$$ clf
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 figw figh])

sc = 1; % slope counter


% during this loop we also compute interannual variables
for i = year(1):year(length(year))
    
    yyyy = str2num(datestr(N_inter,10)); 
    mm = str2num(datestr(N_inter,5));
    dd = str2num(datestr(N_inter,7));

    n = datenum(999, mm, dd); %year 999 for plot
 
    I = find(yyyy==i);

% $$$     figure(1)
% $$$     % find min and its position        
% $$$     plot(mm(I), min(T_inter(:,I)), 'k.');
% $$$ 
% $$$     if i ==  year(1)
% $$$         hold on
% $$$     end
    
    %fit curve
    [P,S]=polyfit(mm(I),  min(T_inter(:,I))', 1); % !!! Xvector is I,% not n(I) !!!
 
    fit = P(2)+mm(I).*P(1);
% $$$     Rplot = plot(mm(I), fit, 'k', 'linewidth', 0.25);
% $$$     set(Rplot, 'color', [0.82 0.82 0.82], 'linewidth', 0.25)
% $$$     
    %keyboard
    % correlation
    [rr, pp] = corrcoef(min(T_inter(:,I)), fit');   
    RR_T(sc) = rr(2);
    PP_T(sc) = pp(2);
    
    %slope
    slope(sc)=P(1); 

    %% - Other interannual variables - %%
         
    % thickness of the CIL
    [YY, JJ] = find(T_inter(:,I)<=1);
    for j=1:length(I)
        %thickness
        thick(j)=length(find(JJ==j));
        
        %heat content
        C = find(JJ==j);
        indi = YY(C);
        DENS = sw_dens(S_inter(indi,I(j)), T_inter(indi,I(j)), Pres(indi)');
             
        if ~isempty(indi)==1
            % - Mean - %
            HH_T(j) = cp*nansum(DENS.*(T_inter(indi, I(j))-freeze_pt))*dz/length(indi);      
            % - Integrated - %
            %HH_T(j) = cp*nansum(DENS.*(T_inter(indi, I(j))-freeze_pt))*dz;              
        else
           HH_T(j) = 0; 
        end
        
    end
    
    %% -- thickness -- %%
    II = find(thick~=0);
    [PP,SS]=polyfit(mm(I(II)),  thick(II)', 1);
    thickness(sc) = PP(1); 
    fit = PP(2)+mm(I(II)).*PP(1); % fit and correlation

    % validation test
% $$$     figure(2)
% $$$     plot(mm(I(II)), thick(II), 'k.')
% $$$     if i ==  year(1)
% $$$         hold on
% $$$     end
% $$$     Rplot = plot(mm(I(II)), fit, 'k', 'linewidth', 0.25);
% $$$     set(Rplot, 'color', [0.82 0.82 0.82], 'linewidth', 0.25)
    
    [rr, pp] = corrcoef(thick(II), fit);   
    RR_D(sc) = rr(2);
    PP_D(sc) = pp(2);
    clear thick

    %% -- heat content -- %%
    Na = find(isnan(HH_T)==1);
    HH_T(Na)=0;
    II = find(HH_T~=0);
    [PP,SS]=polyfit(mm(I(II)),  HH_T(II)', 1);
    HH(sc) = PP(1); 
    fit = PP(2)+mm(I(II)).*PP(1); % fit and correlation
    
    % validation test
% $$$ 
% $$$     figure(3)
% $$$     plot(mm(I(II)), HH_T(II), 'k.')
% $$$     if i ==  year(1)
% $$$         hold on
% $$$     end
% $$$     Rplot = plot(mm(I(II)), fit, 'k', 'linewidth', 0.25);
% $$$     set(Rplot, 'color', [0.82 0.82 0.82], 'linewidth', 0.25)

    [rr, pp] = corrcoef(HH_T(II), fit);   
    RR_H(sc) = rr(2);
    PP_H(sc) = pp(2);
    clear HH_T
    
    sc = sc + 1;

end

% $$$ figure(3)
% $$$ hold off
% $$$ figure(2)
% $$$ hold off
% $$$ 
% $$$ figure(1)
% $$$ hold off
% $$$ ylim([-1.5 2])  
% $$$ xlim(T_LIM);
% $$$ ylab = ylabel('T(^{\circ}C)', 'fontsize', 10, 'VerticalAlignment', 'top');
% $$$ set(gca, 'xticklabel', [], 'fontsize', 10)
% $$$ set(gca, 'xtick', TIK)
% $$$ set(gca, 'fontsize', 10)
% $$$ set(gca, 'XGrid', 'on')
% $$$ 
% $$$ Ypos = get(ylab, 'position');
% $$$ set(ylab, 'position', [Ypos(1)+Ylab_offset Ypos(2) Ypos(3)]);
% $$$ 
% $$$ % -- Output figure costumization -- %
% $$$ Out=get(gca, 'outer'); %get figure box properties
% $$$ Pos=get(gca, 'position');
% $$$ Tig=get(gca, 'tight');
% $$$ 
% $$$ % Increase TightInset for colorbar and Xlabel (forced)
% $$$ Tig = Tig + [Ylab_tightoffset 0 0 0]; 
% $$$ Tig = Tig + [0 offset4 offset5 offset6];
% $$$ % reduce TightInset to reduce white space around figure
% $$$ set(gca, 'Position', Out - Tig * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
% $$$ 
% $$$ set(gcf, 'renderer', 'painters'); % vectorial figure                                  
% $$$ print('-depsc2', 'multiyear_CIL_Tmin_color.eps')                                  
% $$$ %print('-deps2', 'multiyear_CIL_Tmin.eps')
% $$$ 
% $$$ %replace datetick for custumization
% $$$ nn=datenum(999, month, 15);
% $$$ for i=2:length(month)-1
% $$$     mm=datestr(nn(i), 3);
% $$$     day_norm = (nn(i)-T_LIM(1))/no_days;
% $$$     text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
% $$$          'center', 'verticalalignment', 'top', 'fontsize', 10)
% $$$ end
% $$$ print('-depsc2', 'multiyear_CIL_Tmin_xlabel_color.eps')  
%print('-deps2', 'multiyear_CIL_Tmin_xlabel.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Depth Tmin -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 figw figh])
% $$$ 
% $$$ 
% $$$ mc = 1; % mindepth counter
% $$$ sc = 1; % counter
% $$$ 
% $$$ for i = year(1):year(length(year))
% $$$     
% $$$     yyyy = str2num(datestr(N_inter,10)); 
% $$$     mm = str2num(datestr(N_inter,5));
% $$$     dd = str2num(datestr(N_inter,7));
% $$$ 
% $$$     n = datenum(999, mm, dd); %year 999 for plot
% $$$  
% $$$     I = find(yyyy==i);
% $$$ 
% $$$     % find min and its position        
% $$$     [Y, J] = min(T_inter(:,I)); % position and value of Tmin
% $$$ 
% $$$     plot(n(I), J, 'k');
% $$$     
% $$$     plot(n(I), min(T_inter(:,I)), 'k.');
% $$$     
% $$$     if i ==  year(1)
% $$$         hold on
% $$$     end
% $$$ 
% $$$     Zmin(mc:mc+length(Y)-1) = J;
% $$$     Tmin(mc:mc+length(Y)-1) = Y;
% $$$     n_zmin(mc:mc+length(Y)-1) = N_inter(I);
% $$$     mc = mc+length(Y);
% $$$ 
% $$$ end
% $$$ 
% $$$ 
% $$$ hold off
% $$$ 
% $$$ set(gca, 'ydir', 'reverse')
% $$$ 
% $$$ %ylim([-1.5 2])  
% $$$ xlim(T_LIM);
% $$$ ylab = ylabel('{Z_{Tmin}(m)}', 'fontsize', 10, 'VerticalAlignment', 'top');
% $$$ set(gca, 'xticklabel', [], 'fontsize', 10)
% $$$ set(gca, 'xtick', TIK)
% $$$ set(gca, 'fontsize', 10)
% $$$ set(gca, 'XGrid', 'on')
% $$$ 
% $$$ Ypos = get(ylab, 'position');
% $$$ set(ylab, 'position', [Ypos(1)+Ylab_offset Ypos(2) Ypos(3)]);
% $$$ 
% $$$ % -- Output figure costumization -- %
% $$$ Out=get(gca, 'outer'); %get figure box properties
% $$$ Pos=get(gca, 'position');
% $$$ Tig=get(gca, 'tight');
% $$$ 
% $$$ % Increase TightInset for colorbar and Xlabel (forced)
% $$$ Tig = Tig + [Ylab_tightoffset 0 0 0]; 
% $$$ Tig = Tig + [0 offset4 offset5 offset6];
% $$$ % reduce TightInset to reduce white space around figure
% $$$ set(gca, 'Position', Out - Tig * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
% $$$ 
% $$$ set(gcf, 'renderer', 'painters'); % vectorial figure                                  
% $$$ print('-deps2', 'multiyear_CIL_Zmin.eps')
% $$$ 
% $$$ %replace datetick for custumization
% $$$ nn=datenum(999, month, 15);
% $$$ for i=2:length(month)-1
% $$$     mm=datestr(nn(i), 3);
% $$$     day_norm = (nn(i)-T_LIM(1))/no_days;
% $$$     text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
% $$$          'center', 'verticalalignment', 'top', 'fontsize', 10)
% $$$ end

% $$$ print('-deps2', 'multiyear_CIL_Zmin_xlabel.eps')
% $$$ 




% $$$ 
% $$$ 
% $$$ % ----------------------------------------------------------- % 
% $$$ % ----------------- Interannual time series ----------------- %
% $$$ % ----------------------------------------------------------- %
% $$$ 
% Parameters for figure costumization
offset1 = 0.01; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 0.0005; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom
offset5 = 0; % Xtra offset top of figure
offset6 = 0; % Xtra offset right of figure
% $$$ cbar_width = 0.03;
% $$$ clabel_width = 0.08;
% $$$ ti_cbar_frac = 15/16; %reduction of distance bet. colorbar and its
% $$$                       %title
Ylab_offset = -0.8; %offset between yaxis and ylabel
Ylab_tightoffset = 0.09; %TightInset that must be added for title


%%%%%%%%%%%%%%%%%%%%
% -- slope plot -- %
%%%%%%%%%%%%%%%%%%%%


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 figw figh])


plot(year, slope, 'k.')
hold on

% shade area
%y1 = (slope'./slope')*mean(slope) + std(slope);
%y2 = (slope'./slope')*mean(slope) - std(slope);
y1 = (slope'./slope')*0.24+0.04;
y2 = (slope'./slope')*0.24-0.04;
patch([year; flipdim(year,1); year(1)], [y1; flipdim(y2,1); y1(1)], [.8 .8 .8], 'edgecolor', 'none');

plot(year, slope, 'k.')
plot(year, slope, 'k')
plot(year, (slope./slope)*0.24, '--k')

ylim([0.1 0.35])  
xlim([1992 2011]);
%ylab = ylabel('$\dot{T}_{min} (^{\circ}C mo^{-1})$','Interpreter','LaTex', 'fontsize', 10, 'VerticalAlignment', 'top');
ylab = ylabel('T_{min} (^{\circ}C mo^{-1})', 'fontsize', 10, 'VerticalAlignment', 'top');

%set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', 1993:2:2011)
set(gca, 'fontsize', 10)
%set(gca, 'XGrid', 'on')

Ypos = get(ylab, 'position');
set(ylab, 'position', [Ypos(1)+Ylab_offset Ypos(2) Ypos(3)]);

% -- Output figure costumization -- %
Out=get(gca, 'outer'); %get figure box properties
Pos=get(gca, 'position');
Tig=get(gca, 'tight');

% Increase TightInset for colorbar and Xlabel (forced)
Tig = Tig + [Ylab_tightoffset 0 0 0]; 
Tig = Tig + [0 offset4+0.1 offset5 offset6];
% reduce TightInset to reduce white space around figure
set(gca, 'Position', Out - Tig * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
set(gcf, 'renderer', 'painters'); % vectorial figure                                  
                                  
print('-deps2', 'inter_CIL_slope_xlabel.eps')
print('-dpng', '-r300',  'inter_CIL_slope_xlabel.png')

% $$$ set(gca, 'xticklabel', [], 'fontsize', 10)
% $$$ print('-deps2', 'inter_CIL_slope.eps')


keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- CIL thickness / HEAT CONTENT -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 15 4])
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 15 4])

% Compute interannual thickness
for i = 5:2:9 %months on the graph
    
    yyyy = str2num(datestr(N_inter,10)); 
    mm = str2num(datestr(N_inter,5));
    dd = str2num(datestr(N_inter,7));
 
    I = find(mm==i); % vector with all year having i month

    for j = 1:length(I)
        II = find(T_inter(:,I(j))<=T_core);
        d(j) = length(II); % CIL thickness

        DENS = sw_dens(S_inter(II,I(j)), T_inter(II,I(j)), Pres(II)');%density for bin into CIL core
        H(j) = cp*sum(DENS.*(T_inter(II, I(j))-T_core))*dz; % HEAT content
    end
    
    if i == 5
        figure(1)
        plot(yyyy(I), d, '-.k*')        
        hold on
        figure(2)
        plot(yyyy(I), H*1000/1e9, '-.k*')    
        hold on
    end
    if i==7
        figure(1)
        plot(yyyy(I), d, '--ko');
        figure(2)
        plot(yyyy(I), H*1000/1e9, '--ko')
    end
    if i==9
        figure(1)
        plot(yyyy(I), d, ':ks');
        figure(2)
        plot(yyyy(I), H*1000/1e9, ':ks')
    end
    
    
    
end
figure(1)
hold off
%datetick
%legend('May','July','September')

xlim([1992 2010]);
ylab = ylabel('d(m)', 'fontsize', 10, 'VerticalAlignment', 'top');
%set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', 1993:2:2009)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

Ypos = get(ylab, 'position');
set(ylab, 'position', [Ypos(1)+Ylab_offset Ypos(2) Ypos(3)]);

% -- Output figure costumization -- %
Out=get(gca, 'outer'); %get figure box properties
Pos=get(gca, 'position');
Tig=get(gca, 'tight');

% Increase TightInset for colorbar and Xlabel (forced)
Tig = Tig + [Ylab_tightoffset 0 0 0]; 
Tig = Tig + [0 offset4+0.1 offset5 offset6];
% reduce TightInset to reduce white space around figure
set(gca, 'Position', Out - Tig * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

set(gcf, 'renderer', 'painters'); % vectorial figure                                  
print('-deps2', 'inter_CIL_thick_xlabel.eps')

set(gca, 'xticklabel', [], 'fontsize', 10)
print('-deps2', 'inter_CIL_thick.eps')


figure(2)
hold off
%datetick
%legend('May','July','September')

xlim([1992 2010]);
ylim([-0.6 0])
ylab = ylabel('Heat (GJ/m^2)', 'fontsize', 10, 'VerticalAlignment', 'top');
%set(gca, 'xticklabel', [], 'fontsize', 10)
set(gca, 'xtick', 1993:2:2009)
set(gca, 'fontsize', 10)
set(gca, 'XGrid', 'on')

Ypos = get(ylab, 'position');
set(ylab, 'position', [Ypos(1)+Ylab_offset Ypos(2) Ypos(3)]);

% -- Output figure costumization -- %
Out=get(gca, 'outer'); %get figure box properties
Pos=get(gca, 'position');
Tig=get(gca, 'tight');

% Increase TightInset for colorbar and Xlabel (forced)
Tig = Tig + [Ylab_tightoffset 0 0 0]; 
Tig = Tig + [0 offset4+0.1 offset5 offset6];
% reduce TightInset to reduce white space around figure
set(gca, 'Position', Out - Tig * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);



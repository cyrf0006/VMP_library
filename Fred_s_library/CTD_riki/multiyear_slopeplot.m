% just call multiyear_slopeplot
% version from march 2011 with bootstrap on erosion variables
clear

% Parameters to edit

year = [1993:2010]';
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

% Parameters for figure costumization
offset1 = 0.0; % offset between figure and colorbar
offset2 = 0.01; % offset between heigth of colorbar and heigth of figure
offset3 = 0.0005; % offset beteween figure and custom Xlabel 
offset4 = 0.02; % offset between Xlabel and OuterPosition at bottom
offset5 = 0; % Xtra offset top of figure
offset6 = 0; % Xtra offset right of figure
% $$$ cbar_width = 0.03;
% $$$ clabel_width = 0.08;
% $$$ ti_cbar_frac = 15/16; %reduction of distance bet. colorbar and its
% $$$                       %title
%Ylab_offset = -14; %offset between yaxis and ylabel
%Ylab_tightoffset = 0.09; %TightInset that must be added for title

% Param from bin_lenear_itp
depth=1:300;

Tmax=1;
Tmin=-1.5;
VT = Tmin:1:Tmax;

% param for subplot
width = 3;%figures side-by-side
height = round(length(year)/width);
width_each = 5;%cm
height_each = 5;%cm
paperwidth = 15;%cm
paperheight = 20;%cm

% Parameters for figure costumization
offset1 = 0.08; % left of figure 
offset2 = 0.01; % right of figure
offset3 = 0.015; % top of figure
offset4 = -0.004; % bottom of figure (negative because xout has
                  % been floored)
offset5 = 0.015; % between figures dx
offset6 = 0.03; % between figures dy
cbar_width = 0;
cbar_offset = 0; % colorbar offset from figure
                   
%ti_cbar_frac = 1; %reduction of distance bet. colorbar and its
                      
xfigcount=1; %figure counter
yfigcount=1; %id.
xpos_ext=(1-(offset1+(width-1)*offset5+offset2))/width; %position width of each figure
ypos_ext=(1-(offset3+(height-1)*offset6+offset4))/height; %position height of each figure
xout_ext=1/width; %Outerposition width of each figure
yout_ext=1/height;%Outerposition heigth of each figure
% preceeding lines need to be floored
xpos_ext = floor(xpos_ext*100)/100; 
ypos_ext = floor(ypos_ext*100)/100;
xout_ext = floor(xout_ext*100)/100;
yout_ext = floor(yout_ext*100)/100;


% some constant
g = 9.81;
rho_0 = 1.035e3;%kg/m^3
cp = 3.99; %Kj/Kg/K
T_core = 1; %degC, threshold for the CIL core
dz = 1;%m
theta=35;%degrees
Tcil=1;
freeze_pt=-1.8;

Pres = 1:300; %pressure vector


% ----------------------------------------------------------- % 
% ------------------- Compute time serie -------------------- %
% ----------------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%
% -- CIL properties -- %
%%%%%%%%%%%%%%%%%%%%%%%%

% load all casts
T = load('tprofiles.dat');
S = load('sprofiles.dat');
N = load('datprofiles.dat');

yyyy = str2num(datestr(N,10));
mm = str2num(datestr(N,5));
dd = str2num(datestr(N,7));


count = 1; % counter for the number of considered pts
for i = year(1):year(length(year))
    for j=month(1):month(length(month))
        

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

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 paperwidth paperheight])
sc = 1; % slope counter

% during this loop we also compute interannual variables
for i = year(1):year(length(year))
    
    ii = i-year(1)+1;

    yyyy = str2num(datestr(N_inter,10)); 
    mm = str2num(datestr(N_inter,5));
    dd = str2num(datestr(N_inter,7));

    n = datenum(999, mm, dd); %year 999 for plot
 
    I = find(yyyy==i);

    % plot begins
    figure(1)
    subplot(height, width, ii)

    % find min and its position        
    plot(n(I), min(T_inter(:,I)), '--k.');
    hold on
    %fit curve
    [P,S]=polyfit(I,  min(T_inter(:,I))', 1); % !!! Xvector is I,
                                             % not n(I) !!!
    fit = P(2)+I.*P(1);
    Rplot = plot(n(I), fit, 'k', 'linewidth', 1);
    set(Rplot, 'color', 'k', 'linewidth', 1)
    Tmin_slope(sc) = P(1);
    sc = sc+1;
    hold off

    % set output (from bin_linear_itp.m)
    axis([datenum(999, 4,1) datenum(999,11,30) Tmin Tmax])
    %ylim([Tmin Tmax])
    set(gca, 'xticklabel', [], 'fontsize', 10)
    TIK = [datenum(999, month,1)];
    TIK15 = [datenum(999, month,15)];    
    set(gca, 'xtick', TIK)
    set(gca, 'fontsize', 10)      
    set(gca, 'YGrid', 'on')        
    set(gca, 'XGrid', 'on')

    if ii==1 | ii==4 | ii==7 | ii==10 | ii==13 | ii==16
        if ii==7 
            ylabel('T (^{\circ}C)') 
        end
        
    else
        set(gca, 'yticklabel', [], 'fontsize', 10)
    end
    
    if ii==1 | ii==2 | ii==3
        T_LIM = [datenum(999, 4,1) datenum(999,11,30)]; %XLIM for plot
        no_days = T_LIM(2)-T_LIM(1);
        for j=1:length(TIK15)
            mm=datestr(TIK15(j), 4);
            day_norm = (TIK15(j)-T_LIM(1))/no_days;
            text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
                 'center', 'verticalalignment', 'top', 'fontsize', 10)
        end         
    end
    
    
    % Write year
    text(datenum(999, 11, 15), -1.3, sprintf('%d',year(ii)), ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
    
    % ---  Adjust white space between figures! --- %
    % 1) define Outerposition to take all space but no more
    if xfigcount==1 %beginning of a line
        x0 = 0;
    else %not the 1st column
        x0 = (xfigcount-1)*(xout_ext);
    end
    
    y0 = 1 - yout_ext*yfigcount;
    
    Out = [x0 y0 xout_ext yout_ext];
    
    % 2) define position for minimizing white space
    if xfigcount==1 %beginning of a line
        x0 = offset1;
    else %not the 1st column
        x0 = offset1+(xfigcount-1)*(xpos_ext+offset5);
    end
    
    if yfigcount==1  %1st line
        y0 = 1 - offset3 - ypos_ext ;
    else  %not the 1st line
        y0 = 1 - (offset3 + (yfigcount-1)*offset6 + yfigcount*ypos_ext);
    end

    Pos = [x0 y0  xpos_ext  ypos_ext];
    set(gca, 'XGrid', 'on')

    set(gca, 'outerposition', Out) % set outerposition
    set(gca, 'position', Pos) % set position
    % ------------- end adjust space -------------------- %

        % increment subplot counter (for either T or S)
    if xfigcount<width
        xfigcount = xfigcount+1;
    else
       xfigcount = 1; 
       yfigcount = yfigcount+1;
    end

end

set(gcf, 'renderer', 'painters'); % vectorial figure                                  
print('-depsc2', 'multiyear_CIL_Tmin_sub.eps')                                  





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- CIL thickness and heat content -- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reset counters
xfigcount=1; %figure counter
yfigcount=1; %id.
xpos_ext=(1-(offset1+(width-1)*offset5+offset2))/width; %position width of each figure
ypos_ext=(1-(offset3+(height-1)*offset6+offset4))/height; %position height of each figure
xout_ext=1/width; %Outerposition width of each figure
yout_ext=1/height;%Outerposition heigth of each figure
% preceeding lines need to be floored
xpos_ext = floor(xpos_ext*100)/100; 
ypos_ext = floor(ypos_ext*100)/100;
xout_ext = floor(xout_ext*100)/100;
yout_ext = floor(yout_ext*100)/100;



figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 paperwidth ...
                    paperheight])

figure(3)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 paperwidth ...
                    paperheight])

sc = 1; % slope counter

% during this loop we also compute interannual variables
for i = year(1):year(length(year))
    
    ii = i-year(1)+1;

    yyyy = str2num(datestr(N_inter,10)); 
    mm = str2num(datestr(N_inter,5));
    dd = str2num(datestr(N_inter,7));

    n = datenum(999, mm, dd); %year 999 for plot
 
    I = find(yyyy==i);

    clear thickness heatcontent
    for j = 1:length(I)
    
        CIL =  find(T_inter(:,I(j))<Tcil);
        if ~isempty(CIL) == 1
            DENS = sw_dens(S_inter(CIL,I(j)), T_inter(CIL,I(j)), depth(CIL)'); 
            heatcontent(j) =  cp*nansum(DENS.*(T_inter(CIL, I(j))-freeze_pt))*dz/length(CIL)/1000;%MJ
            thickness(j) = length(CIL);  
        else
            heatcontent(j) = NaN;  
            thickness(j) = 0;
        end
    
    end
    thickness = thickness';
    heatcontent = heatcontent';
    
    
    figure(2)
    subplot(height, width, ii)
    % find min and its position        
    plot(n(I), thickness, '--k.');
    hold on
    %fit curve
    [P,S]=polyfit(I,  thickness, 1); % !!! Xvector is I,
                                             % not n(I) !!!
    fit = P(2)+I.*P(1);
    Rplot = plot(n(I), fit, 'k', 'linewidth', 1);
    set(Rplot, 'color', 'k', 'linewidth', 1)
    thick_slope(sc) = P(1);
    hold off
    
    figure(3)
    subplot(height, width, ii)
    % find min and its position        
    plot(n(I), heatcontent, '--k.');
    hold on
    %fit curve
    II = find(~isnan(heatcontent)==1);
    [P,S]=polyfit(I(II),  heatcontent(II), 1); % !!! Xvector is I,
                                             % not n(I) !!!
    fit = P(2)+I(II).*P(1);
    Rplot = plot(n(I(II)), fit, 'k', 'linewidth', 1);
    set(Rplot, 'color', 'k', 'linewidth', 1)
    heatc_slope(sc) = P(1);
    hold off

    figure(2)
    % set output (from bin_linear_itp.m)
    axis([datenum(999, 4,1) datenum(999,11,30) 0 100])
    %ylim([Tmin Tmax])
    set(gca, 'xticklabel', [], 'fontsize', 10)
    TIK = [datenum(999, month,1)];
    TIK15 = [datenum(999, month,15)];    
    set(gca, 'xtick', TIK)
    set(gca, 'fontsize', 10)      
    set(gca, 'YGrid', 'on')        
    set(gca, 'XGrid', 'on')

    if ii==1 | ii==4 | ii==7 | ii==10 | ii==13 | ii==16
        if ii==7 
            ylabel('T (^{\circ}C)') 
        end
        
    else
        set(gca, 'yticklabel', [], 'fontsize', 10)
    end
    
    if ii==1 | ii==2 | ii==3
        T_LIM = [datenum(999, 4,1) datenum(999,11,30)]; %XLIM for plot
        no_days = T_LIM(2)-T_LIM(1);
        for j=1:length(TIK15)
            mm=datestr(TIK15(j), 4);
            day_norm = (TIK15(j)-T_LIM(1))/no_days;
            text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
                 'center', 'verticalalignment', 'top', 'fontsize', 10)
        end         
    end
    
    
    % Write year
    text(datenum(999, 11, 15), -1.3, sprintf('%d',year(ii)), ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
    
    % ---  Adjust white space between figures! --- %
    % 1) define Outerposition to take all space but no more
    if xfigcount==1 %beginning of a line
        x0 = 0;
    else %not the 1st column
        x0 = (xfigcount-1)*(xout_ext);
    end
    
    y0 = 1 - yout_ext*yfigcount;
    
    Out = [x0 y0 xout_ext yout_ext];
    
    % 2) define position for minimizing white space
    if xfigcount==1 %beginning of a line
        x0 = offset1;
    else %not the 1st column
        x0 = offset1+(xfigcount-1)*(xpos_ext+offset5);
    end
    
    if yfigcount==1  %1st line
        y0 = 1 - offset3 - ypos_ext ;
    else  %not the 1st line
        y0 = 1 - (offset3 + (yfigcount-1)*offset6 + yfigcount*ypos_ext);
    end

    Pos = [x0 y0  xpos_ext  ypos_ext];
    set(gca, 'XGrid', 'on')

    set(gca, 'outerposition', Out) % set outerposition
    set(gca, 'position', Pos) % set position
    % ------------- end adjust space -------------------- %
    
    figure(3)
    % set output (from bin_linear_itp.m)
    axis([datenum(999, 4,1) datenum(999,11,30) 6 12])
    %ylim([Tmin Tmax])
    set(gca, 'xticklabel', [], 'fontsize', 10)
    TIK = [datenum(999, month,1)];
    TIK15 = [datenum(999, month,15)];    
    set(gca, 'xtick', TIK)
    set(gca, 'fontsize', 10)      
    set(gca, 'YGrid', 'on')        
    set(gca, 'XGrid', 'on')

    if ii==1 | ii==4 | ii==7 | ii==10 | ii==13 | ii==16
        if ii==7 
            ylabel('T (^{\circ}C)') 
        end
        
    else
        set(gca, 'yticklabel', [], 'fontsize', 10)
    end
    
    if ii==1 | ii==2 | ii==3
        T_LIM = [datenum(999, 4,1) datenum(999,11,30)]; %XLIM for plot
        no_days = T_LIM(2)-T_LIM(1);
        for j=1:length(TIK15)
            mm=datestr(TIK15(j), 4);
            day_norm = (TIK15(j)-T_LIM(1))/no_days;
            text(day_norm,-offset3,mm,'units','normalized', 'HorizontalAlignment', ...
                 'center', 'verticalalignment', 'top', 'fontsize', 10)
        end         
    end
    
    
    % Write year
    text(datenum(999, 11, 15), -1.3, sprintf('%d',year(ii)), ...
         'verticalalignment', 'bottom', 'horizontalalignment', 'right', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);
    
    % ---  Adjust white space between figures! --- %
    % 1) define Outerposition to take all space but no more
    if xfigcount==1 %beginning of a line
        x0 = 0;
    else %not the 1st column
        x0 = (xfigcount-1)*(xout_ext);
    end
    
    y0 = 1 - yout_ext*yfigcount;
    
    Out = [x0 y0 xout_ext yout_ext];
    
    % 2) define position for minimizing white space
    if xfigcount==1 %beginning of a line
        x0 = offset1;
    else %not the 1st column
        x0 = offset1+(xfigcount-1)*(xpos_ext+offset5);
    end
    
    if yfigcount==1  %1st line
        y0 = 1 - offset3 - ypos_ext ;
    else  %not the 1st line
        y0 = 1 - (offset3 + (yfigcount-1)*offset6 + yfigcount*ypos_ext);
    end

    Pos = [x0 y0  xpos_ext  ypos_ext];
    set(gca, 'XGrid', 'on')

    set(gca, 'outerposition', Out) % set outerposition
    set(gca, 'position', Pos) % set position
    % ------------- end adjust space -------------------- %

    
    % increment subplot counter (for either thickness and heatcontent)
    if xfigcount<width
        xfigcount = xfigcount+1;
    else
       xfigcount = 1; 
       yfigcount = yfigcount+1;
    end 
    
    
    %slope counter
    sc = sc+1;

end


% $$$ figure(1)
% $$$ set(gcf, 'renderer', 'painters'); % vectorial figure                            
% $$$ print('-depsc2', 'multiyear_CIL_Tmin_sub.eps')   
% $$$ 

% -- bootstrap on erosion variables -- %
disp('bootstrap...')
nboot = 1000;
N = length(Tmin_slope)

for b = 1:nboot
    % slope
    r = rand(N,1);
    r = ceil(r*N/1);
    tmin_slope_boot_b(b) = nanmean(Tmin_slope(r),2);
    thick_slope_boot_b(b) = nanmean(thick_slope(r),2);
    heatc_slope_boot_b(b) = nanmean(heatc_slope(r),2);
end

% Compute average
tmin_slope_boot_dot = nanmean(tmin_slope_boot_b);
thick_slope_boot_dot = nanmean(thick_slope_boot_b);
heatc_slope_boot_dot = nanmean(heatc_slope_boot_b);

% Compute standard error
tmin_slope_error = sqrt(sum((diff([tmin_slope_boot_b tmin_slope_boot_dot],1, 2)).^2, 2)./(nboot-1));
thick_slope_error = sqrt(sum((diff([thick_slope_boot_b thick_slope_boot_dot],1, 2)).^2, 2)./(nboot-1));
heatc_slope_error = sqrt(sum((diff([heatc_slope_boot_b heatc_slope_boot_dot],1, 2)).^2, 2)./(nboot-1));

% Compute 95% confidence interval
tmin_sort = sort(tmin_slope_boot_b, 2);
thick_sort  = sort(thick_slope_boot_b, 2);    
heatc_sort = sort(heatc_slope_boot_b,2);

CI_2p5 = round(2.5/100*nboot)
CI_97p5 = round(97.5/100*nboot)

tmin_2p5 = tmin_sort(:,CI_2p5);
thick_2p5 = thick_sort(:,CI_2p5);
heatc_2p5 = heatc_sort(:,CI_2p5);
tmin_97p5 = tmin_sort(:,CI_97p5);
thick_97p5 = thick_sort(:,CI_97p5);
heatc_97p5 = heatc_sort(:,CI_97p5);

disp('Tmin...')
[tmin_slope_boot_dot tmin_2p5 tmin_97p5]
disp(' ')
disp('thickness...')
[thick_slope_boot_dot thick_2p5 thick_97p5]
disp(' ')
disp('heat content...')
[heatc_slope_boot_dot heatc_2p5 heatc_97p5]



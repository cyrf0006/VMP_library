clear

%load BR_colormap
%load BRcolormap3 % noblue under 1degC
load('/home/cyrf0006/PhD/CTD_IML4/BIN/renamed/BRcolormap3')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- parameters to edit --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ctd_files = 'ctd_files';
month = 4:11;
yearnum=1991:2008; %this could be edited
                   %yearstr=datestr(dat, 10);
depth=1:300;

Tmax=6;
Tmin=1.5;
Smax=34;
Smin=25;
VT = Tmin:0.1:Tmax;
VS = Smin:1:Smax;

% param for subplot
width = 3;%figures side-by-side
height = round(length(yearnum)/width);
width_each = 5;%cm
height_each = 5;%cm
paperwidth = 14;%cm                
paperheight = 17;%cm
%paperheight = 15;%cm


% Parameters for figure costumization
offset1 = 0.09; % left of figure 
offset2 = 0.005; % right of figure
                 
%offset3 = 0.046; % top of figure
offset3 = 0.055;
offset4 = 0.09 ;                
%offset4 = -0.025; % bottom of figure (negative because xout has
                  % been floored)
offset5 = 0.015; % between figures dx
offset6 = 0.04; % between figures dy

dateoffset = 1.31;
cbar_width = 0.03;
%cbar_offset = 0.01; % colorbar offset from figure
cbar_offset = 0.005;
ti_cbar_frac = 1; %reduction of distance bet. colorbar and its

offset10 = 0.02; %offset between Pos and Out in y

xfigcount=1; %figure counter
yfigcount=1; %id.
xpos_ext=(1-(offset1+(width-1)*offset5+offset2))/width; %position width of each figure
ypos_ext=(1-(offset3+(height-1)*offset6+offset4))/height; %position height of each figure
xout_ext=1/width;

 %Outerposition width of each figure
yout_ext=1/height;%Outerposition heigth of each figure
% preceeding lines need to be floored
xpos_ext = floor(xpos_ext*100)/100; 
ypos_ext = floor(ypos_ext*100)/100;
xout_ext = floor(xout_ext*100)/100;
yout_ext = floor(yout_ext*100)/100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- load data (raw profile, listed in ctd_files) --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ n = load('datprofiles.dat');
% $$$ Tmat = load('tprofiles.dat');
% $$$ Smat = load('sprofiles.dat');

[Tmat, Smat, n] = ctd_matrix_baltic(ctd_files);


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paperwidth paperheight])

% $$$ figure(2)
% $$$ clf
% $$$ set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 paperwidth paperheight])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- loop on year to plot ---- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loop on years')
count=0;
for i = 1:length(yearnum)
    
    disp(sprintf('  %d', yearnum(i)))
    % new daily axis
    time_axis = datenum(yearnum(i), min(month), 1): datenum(yearnum(i), ...
                                                      max(month), 30); ...
    
    depth_axis = depth;  
    
    % select profiles in considered year and in month range
    I = find(str2num(datestr(n,10))==yearnum(i) & ...
                     str2num(datestr(n,5))>=min(month) & ...
                     str2num(datestr(n,5))<=max(month));
      
    Tyear = Tmat(:,I);
    Syear = Smat(:,I);
    n_year = n(I);
  
    %remove profiles when more than one / day
    [B, I] = unique(round(n_year));
    n_year = n_year(I);
    Tyear = Tyear(:,I);
    Syear = Syear(:,I);   
  
    if size(Tyear, 2) < 2
        continue
    else
        count=count+1;
    end
    
    
    tprofile_tick = n_year; % to identify profile on contourplot
    sprofile_tick = n_year;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- loop on bin ---- %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    % Here we consider all good data for any year and store them in
    % vectors
    for j=depth       
        TJ=find(isnan(Tyear(j,:))==0); %ignore nans
        SJ=find(isnan(Syear(j,:))==0); %ignore nans
        if length(TJ)>1
            T_itp(j,:) = interp1(n_year(TJ), Tyear(j, TJ), time_axis);
        else
            T_itp(j,1:length(time_axis))=NaN;
        end
        if length(SJ)>1
            S_itp(j,:) = interp1(n_year(SJ), Syear(j, SJ), time_axis);
        else
            S_itp(j,:)=NaN;
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% -- Temperature plot -- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    figure(1)
    
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
    
    Pos = [x0 y0+offset10  xpos_ext  ypos_ext]; % it appears that
                                                % offset10 is
                                                % important!
    
    % ------------- end adjust space -------------------- %
        
    % Subplot calling
    subplot(height, width, i, 'outerposition', Out, 'position', Pos);

    % Contour plot
    contourf(time_axis, depth_axis, T_itp, VT, 'linestyle', 'none');
    
    % Xtra contour around CIL
    hold on
    contour(time_axis, depth_axis, T_itp, [1 1 1], 'edgecolor', [0 0 0], 'LineStyle', '-','linewidth', 0.25 )
    hold off
        
    colormap(mycolormap)
    caxis([Tmin, Tmax]);
    axis([datenum(yearnum(i), 4,1) datenum(yearnum(i),11,30) 0 300])
    set(gca, 'ydir', 'reverse')
    set(gca, 'xticklabel', [], 'fontsize', 10)
    TIK = [datenum(yearnum(i), month,1)];
    TIK15 = [datenum(yearnum(i), month,15)];    
    set(gca, 'xtick', TIK)
    set(gca, 'fontsize', 10)
    set(gca, 'XGrid', 'on')
    set(gca, 'YGrid', 'on')
    
    if count==1 | count==4 | count==7 | count==10 | count==13 | count==16
        if count==7 
            ylabel('P(dbar)') 
        end
        
    else
        set(gca, 'yticklabel', [], 'fontsize', 10)
    end
    
    if count==1 | count==2 | count==3
        T_LIM = [datenum(yearnum(i), 4,1) datenum(yearnum(i),11,30)]; %XLIM for plot
        no_days = T_LIM(2)-T_LIM(1);
        for j=1:length(TIK15)
            mm=datestr(TIK15(j), 4);
            day_norm = (TIK15(j)-T_LIM(1))/no_days;
            text(day_norm,dateoffset,mm,'units','normalized', 'HorizontalAlignment', 'center', 'verticalalignment', 'top', 'fontsize',10)
            
        end         
    end
    
    
    % Write year
    text(datenum(yearnum(i), 4, 4), 298, sprintf('%d',yearnum(i)), ...
         'verticalalignment', 'bottom', 'fontsize', 10, 'fontweight', ...
         'bold','BackgroundColor',[1 1 1]);

   
    % Tick where there are profile
    for j=1:length(tprofile_tick)
        %text(tprofile_tick(j), 0, {'\bullet'}, 'verticalalignment', ...
        %     'baseline', 'horizontalalignment', 'center');
        text(tprofile_tick(j), 0, 'l', 'verticalalignment', ...
             'baseline', 'horizontalalignment', 'center', 'fontsize', ...
             6, 'fontweight', 'bold');
        
    end
    
    % set colorbar and its title
    if i==length(yearnum)-1;
        c = colorbar('location', 'southoutside', 'fontsize', 10);
        set(c, 'position', [Pos(1)+cbar_offset Pos(2)-cbar_offset-cbar_width Pos(3)-2*cbar_offset cbar_width*.75]);
        set(c, 'tickdir', 'out')
        %set(c, 'xtick', [-2:5])
        ti = ylabel(c,'{T(^{\circ}C)}', 'FontSize', 10);
        ti_pos = get(ti, 'position');
        
        
        set(ti, 'rotation', 0)
        set(ti, 'position', [4 ti_pos(2)-2 ti_pos(3)]); 
                
    end
        
% $$$     set(gca, 'outerposition', Out) % set outerposition
% $$$     set(gca, 'position', Pos) % set position 
% $$$     
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     %% --- Salinity plot --- %
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     figure(2)
% $$$     set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 paperwidth paperheight])
% $$$     subplot(height, width, i)
% $$$     contourf(time_axis, depth_axis, S_itp, VS);
% $$$     colormap(mycolormap)
% $$$     caxis([Smin Smax]);
% $$$     axis([datenum(yearnum(i), 4,1) datenum(yearnum(i),11,30) 0 300])
% $$$     set(gca, 'ydir', 'reverse')
% $$$     set(gca, 'xticklabel', [], 'fontsize', 10)
% $$$     TIK = [datenum(yearnum(i), month,1)];
% $$$     TIK15 = [datenum(yearnum(i), month,15)];    
% $$$     set(gca, 'xtick', TIK)
% $$$     set(gca, 'fontsize', 10)
% $$$     set(gca, 'XGrid', 'on')
% $$$     set(gca, 'YGrid', 'on')
% $$$     
% $$$     if i==1 | i==4 | i==7 | i==10 | i==13 | i==16
% $$$         if i==7 
% $$$             ylabel('depth (m)') 
% $$$         end
% $$$         
% $$$     else
% $$$         set(gca, 'yticklabel', [], 'fontsize', 10)
% $$$     end
% $$$     
% $$$     if i==1 | i==2 | i==3
% $$$         T_LIM = [datenum(yearnum(i), 4,1) datenum(yearnum(i),11,30)]; %XLIM for plot
% $$$         no_days = T_LIM(2)-T_LIM(1);
% $$$         for j=1:length(TIK15)
% $$$             mm=datestr(TIK15(j), 4);
% $$$             day_norm = (TIK15(j)-T_LIM(1))/no_days;
% $$$             text(day_norm,dateoffset,mm,'units','normalized', 'HorizontalAlignment', ...
% $$$                  'center', 'verticalalignment', 'top', 'fontsize', 10)
% $$$         end         
% $$$     end
% $$$     
% $$$     
% $$$     % Write year
% $$$     text(datenum(yearnum(i), 4, 4), 298, sprintf('%d',yearnum(i)), ...
% $$$          'verticalalignment', 'bottom', 'fontsize', 10, 'fontweight', ...
% $$$          'bold','BackgroundColor',[1 1 1]);
% $$$ 
% $$$     % Tick where there are profile
% $$$     for j=1:length(sprofile_tick)
% $$$         text(sprofile_tick(j), 0, {'\bullet'}, 'verticalalignment', 'baseline');
% $$$     end
% $$$     
% $$$     % ---  Adjust white space between figures! --- %
% $$$     % 1) define Outerposition to take all space but no more
% $$$     if xfigcount==1 %beginning of a line
% $$$         x0 = 0;
% $$$     else %not the 1st column
% $$$         x0 = (xfigcount-1)*(xout_ext);
% $$$     end
% $$$     
% $$$     y0 = 1 - yout_ext*yfigcount;
% $$$     
% $$$     Out = [x0 y0 xout_ext yout_ext];
% $$$     
% $$$     % 2) define position for minimizing white space
% $$$     if xfigcount==1 %beginning of a line
% $$$         x0 = offset1;
% $$$     else %not the 1st column
% $$$         x0 = offset1+(xfigcount-1)*(xpos_ext+offset5);
% $$$     end
% $$$     
% $$$     if yfigcount==1  %1st line
% $$$         y0 = 1 - offset3 - ypos_ext ;
% $$$     else  %not the 1st line
% $$$         y0 = 1 - (offset3 + (yfigcount-1)*offset6 + yfigcount*ypos_ext);
% $$$     end
% $$$ 
% $$$     Pos = [x0 y0  xpos_ext  ypos_ext];
% $$$ 
% $$$     set(gca, 'outerposition', Out) % set outerposition
% $$$     set(gca, 'position', Pos) % set position
% $$$     % ------------- end adjust space -------------------- %
% $$$ 
% $$$     % set colorbar and its title
% $$$     if i==length(yearnum);
% $$$         c = colorbar('FontSize', 10, 'position', [Pos(1)+Pos(3)+cbar_offset ...
% $$$                             Pos(2)+cbar_offset cbar_width Pos(4)-2*cbar_offset]);
% $$$         ti = ylabel(c,'psu', 'FontSize', 10);
% $$$         ti_pos = get(ti, 'position');
% $$$         set(ti, 'position', [ti_pos(1)*ti_cbar_frac ti_pos(2) ti_pos(3)]); 
% $$$         set(gca, 'position', Pos) % reajust position  
% $$$         
% $$$     end
    
    %%%%%%%%%%%%% End plotting %%%%%%%%%%%%%%%%%%%
    
    
    % increment subplot counter (for either T or S)
    if xfigcount<width
        xfigcount = xfigcount+1;
    else
       xfigcount = 1; 
       yfigcount = yfigcount+1;
    end
    
    
% $$$     %save interpolated matrix...
% $$$     Tfname = sprintf('T_itp%d.dat', yearnum(i));
% $$$     Sfname = sprintf('S_itp%d.dat', yearnum(i));
% $$$     Nfname = sprintf('N_itp%d.dat', yearnum(i));
% $$$  
% $$$     dlmwrite(Tfname, T_itp,'delimiter',' ','precision',6)
% $$$     dlmwrite(Sfname, S_itp,'delimiter',' ','precision',6)
% $$$     dlmwrite(Nfname, time_axis,'delimiter',' ','precision',6)

    
    
    %clear matrix
    clear T_itp S_itp %n_year time_axis depth_axis


end
 
disp('Save figures')
    % Save figures
    figure(1)
    set(gcf, 'renderer', 'painters'); % vectorial figure
                                      
    print('-depsc2', 'T_interannual.eps');
    print('-dpng', '-r300', 'T_interannual.png');
    
% $$$     figure(2)
% $$$     set(gcf, 'renderer', 'painters'); % vectorial figure
% $$$     print('-depsc2', 'S_interannual.eps');  
     
function plot_dissip(files)
% plot_dissip(files)
%
% Where files contains the name of the dissip_matrix, one or many ex:
% files = ['dissip_matrix_21-07'; 'dissip_matrix_23-07'; 'dissip_matrix_27-07'; ... ];
% ";" are important!!
%
% This function plot the dissipation for each overturn in a
% depth/time-vs-HighTide figure
% The script acts in 2 parts: building a M.mat matrix and "imagesc-it". 
%
% If you are satistied with a M.mat matrix already, you can only execute
% the second part by giving an empty "files" as parameter, i.e.:
% plot_dissip([])
%
% files =['dissip_matrix_07-21.mat';'dissip_matrix_07-23.mat';'dissip_matrix_07-27.mat';'dissip_matrix_09-09.mat';'dissip_matrix_09-10.mat';'dissip_matrix_09-16.mat']; 
%
% Author: Frederic Cyr - 2009/01/12
%
%
%% -------------------------------------------------------------- %%




if ~isempty(files)
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %     first part     %
    % build the M matrix %
    %%%%%%%%%%%%%%%%%%%%%%
    
    dim= size(files);
    no_files = dim(1);

    if no_files > 1 % if more than one matrix to plot

        build_M(files(1, :), 0)

        for i=2:no_files
            build_M(files(i, :), 1)
        end

    else

        build_M(files(1, :), 0) % just one matrix (i.e. plot only one day)

    end



    %%%%%%%%%%%%%%%%%%%%%
    %    second part    %
    % plot the M matrix %
    %%%%%%%%%%%%%%%%%%%%%

    clear
    load M

    depth = 0:0.1:180;
    depth = single(depth); % 2 digits precision for the following find
    tide = -6:0.1:6;
    tide = single(tide); % 2 digits precision for the following find

    % dissipation matrix
     %CLIM = [-10 -3];
     %imagesc(tide, depth, log10(mat), CLIM)
     imagesc(tide, depth, log10(mat))
     % tide representation
     y = -60*cos(pi.*tide./12)+180;
    hold on
    plot(tide, y, 'r')
    hold off
    %title('Dissipation turbulente pour chaque "overturn" (log10(epsilon))', 'FontSize', 15);
    xlabel('temps par rapport a la maree haute (h)', 'FontSize', 15)
    ylabel('Profondeur (m)', 'FontSize', 15)
    legend('maree', 'FontSize', 15)
    c = colorbar;
    ylabel(c,'{log \epsilon (W/Kg)}', 'FontSize', 15)
    set(gca, 'FontSize', 15)

else  %Only second part because M.mat exists
    
    %%%%%%%%%%%%%%%%%%%%%
    %    second part    %
    % plot the M matrix %
    %%%%%%%%%%%%%%%%%%%%%

    clear
    load M

    depth = 0:0.1:180;
    depth = single(depth); % 2 digits precision for the following find
    tide = -6:0.1:6;
    tide = single(tide); % 2 digits precision for the following find

    % dissipation matrix
     %CLIM = [-10 -3];
     %imagesc(tide, depth, log10(mat), CLIM)
     imagesc(tide, depth, log10(mat))
     % tide representation
     y = -60*cos(pi.*tide./12)+180;
    hold on
    plot(tide, y, 'r')
    hold off
    
    title('Dissipation turbulente pour chaque "overturn"', 'FontSize', 15);
    xlabel('temps p/r maree haute (h)', 'FontSize', 15)
    ylabel('profondeur (m)', 'FontSize', 15)
    legend('maree', 'FontSize', 15, 'Location','SouthWest')
    c = colorbar
    ylabel(c,'log(\epsilon) (W/Kg)', 'FontSize', 15)
    set(gca, 'FontSize', 15)  
    
    
end


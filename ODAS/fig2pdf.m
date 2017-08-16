%% fig2pdf
% Convert figures into PDF files with better quality then exporting
% from Matlab.
%%
% <latex>\index{Type B!fig2pdf}</latex>
%
%%% Syntax
%   fig2pdf( 'figure', 'fileName', [dimensions], legendPosition )
%
% * [figure]   Either a string representing a .fig file or a figure object.
% * [fileName] Output eps/pdf file names.  Only used when 'figure' is a 
%              figure object.
% * [dimensions] Width and height meadured in inches.  Optional, omit or leave
%              empty to use the default value of [8,6].
% * [legendPosition] Position of the legend.  An empty ('') string uses the 
%              position defined by the .fig file.  Default value = 
%              'NorthEastOutside'
%
%%% Description
% Convert a figure object/file to an appropriately scaled PDF 
% file.  Avoids use of the Matlab PDF export function as it does not generate
% attractive plots.  Instead, this function exports to an EPS file then 
% converts the EPS file into the final PDF file.  The PDF file maintains the 
% same bounding box as the exported EPS file.
%
% Axes set to font size 10 and title set to font size 12.  To make the 
% fonts larger/smaller - adjust the dimensions accordingly.  Smaller 
% dimensions make the fonts appear larger.
%
% Legend is positioned in the 'OutsideNorthEast' corner by default.
%
% The resulting graph is saved as a vector so you don't have to worry 
% about the dimensions not being correct.  Image will be attractive even
% when scaled.  Plotting is performed with the "-painters" option.
%
% Resulting PDF files have the same name as figName.  Old files will be 
% overwritten. Figures can have multiple sub-plots so long as they are
% arranged vertically.
%
% Note: LaTeX should be installed on the computer.  This function requires 
% access to the "ps2pdf" program included with most LaTeX distributions.
%
%%% Examples:
%
%    >> figure
%    >> x = 0:0.1:8*pi;
%    >> plot( [sin(x); cos(x)] )
%    >> fig2pdf( gcf, 'my_plot' )
%
% Generating vector images (PDF) from the current figure.  The final output
% will look the same regardless of the platform used to generate the plot.
%
%    >> fig2pdf( 'some_figure_file' )
%
%    >> fig2pdf( 'some_figure_file', [8,4], 'SouthEastInside' )
%
%    >> fig2pdf( 'some_figure_file', [], 'SouthEastInside' )
%
% Generate PDF versions of the figure named 'some_figure_file.fig'.  
% Use of '.fig' files allows this function to be used in scripts.  Because the
% function produces the same images regardless of platform, the scripts simplify
% working with multiple platforms.
%
% The second function call modifies the image size and legend position while 
% the third fuction call modifies the legend position but uses the default [8,6]
% image size.

% *Version History:*
% 2012-08-13 (WID) Enabled function to work with figure objects.  Added
%                  documentation to allow integration into the ODAS library.
% 2012-09-15 (WID) Added support for subplots.  Updated documentation.
% 2012-11-05 (WID) Updated documentation - changed name to fig2pdf.
% 2012-11-07 (WID) Subplots are now all the same width - fixed for quick_bench.
% 2012-11-13 (WID) Extended operation to allow for horizontal and verticle
%                  subplots.
% 2012-12-14 (WID) Major changes to improve performance with different graphs
% 2012-12-14 (WID) Will now save EPS files when ps2pdf is not found

% REFERENCE:  Acceptable Legend Locations (from Matlab documentation)
% North - Inside plot box near top
% South - Inside bottom
% East - Inside right
% West - Inside left
% NorthEast - Inside top right (default for 2-D plots)
% NorthWest - Inside top left
% SouthEast - Inside bottom right
% SouthWest - Inside bottom left
% NorthOutside - Outside plot box near top
% SouthOutside - Outside bottom
% EastOutside - Outside right
% WestOutside - Outside left
% NorthEastOutside - Outside top right (default for 3-D plots)
% NorthWestOutside - Outside top left
% SouthEastOutside - Outside bottom right
% SouthWestOutside - Outside bottom left
% Best - Least conflict with data in plot
% BestOutside - Least unused space outside plot

% HINTS:
% The following returns the current active figure:
% fig = get(0,'CurrentFigure');
% fig = gcf;

function fig2pdf(fig_in, varargin)

% User adjustable default values.
% Dimentions, increase dimensions to reduce font size.
DEFAULT_dimensions = [8,6];
% Axis font size must be at least 10 to fit inside the label.
DEFAULT_axisFontSize = 10;
% Position of legend.
DEFAULT_legendPosition = 'NorthEastOutside';
% A small spacer - a couple of points are all that is required.
DEFAULT_spacer = 0.01;
% *************************************************************


% Set default values for variable arguments...
legendPosition = DEFAULT_legendPosition;
dimensions = DEFAULT_dimensions;

% We have two options.  First, a figure object is sent in so the second
% argument will be the destination eps/pdf file name.  Second, it is a
% file so the second argument is not used.

input_ptr = 1;
use_file = 1;
figName = fig_in;

if ~isa(fig_in, 'char'),
  
  use_file = 0;
  fig = fig_in;

  if nargin < 2,
    fig = get(0, 'CurrentFigure');
    figName = fig_in;
  else
    figName = varargin{input_ptr};
    input_ptr = input_ptr + 1;
  end
end

[P,N,E]=fileparts(figName);

figName = fullfile(P, [N '.fig']);
epsName = fullfile(P, [N '.eps']);
pdfName = fullfile(P, [N '.pdf']);

if length(varargin) >= input_ptr + 1,
  legendPosition = varargin{input_ptr + 1};
end
if length(varargin) >= input_ptr,
  dimensions = varargin{input_ptr}; 
  if isempty(dimensions), dimensions = DEFAULT_dimensions; end
end

position = [0, 0, dimensions(1), dimensions(2)];

if use_file, fig = open(figName); end


% Make a backup copy of the figure
%fig = CloneFig(fig);

% Set the default font sizes.
%set( fig, 'DefaulttextFontSize', DEFAULT_axisFontSize );

set( fig, 'PaperUnits', 'inches' );
set( fig, 'PaperSize', dimensions );
set( fig, 'Units', 'inches' );
set( fig, 'Position', [0.01, 0.01, dimensions(1)-0.01, dimensions(2)-0.01] );
set( fig, 'renderer', 'painters' );

for child = get(fig,'Children')',
    set(child, 'FontSize', DEFAULT_axisFontSize);
    set(get(child,'XLabel'), 'FontSize', DEFAULT_axisFontSize);
    set(get(child,'YLabel'), 'FontSize', DEFAULT_axisFontSize);

    % Set the title to a slightly larger font size.
    set(get(child, 'Title'), 'FontSize', DEFAULT_axisFontSize * 1.2);
end

% Set the legend positions
position_legends;

% Now determine if the subplots are arranged horizontally or vertically and
% position accordingly.  To determin this, find the center of each plot and
% measure the differences between the x and y positions of the center points.

subplots = findobj(fig,'Type','axes','-not','Tag','legend')';
positions = [];
for child = subplots,
    positions(end+1,:) = get(child, 'position');
end

centers = zeros(length(subplots), 2);
for ii = 1:length(subplots),    
    centers(ii,1) = positions(ii,1) + positions(ii,3) / 2;
    centers(ii,2) = positions(ii,2) + positions(ii,4) / 2;
end

difference = std(centers);

if length(difference) < 2,
    % There is only one plot
    arrange_verticle;
elseif difference(1) > difference(2),
    % The plots are arranged horizontally
    arrange_horizontal;
else
    % The plots are arranged vertically
    arrange_verticle;
end

print('-depsc2', epsName)
%close(fig);

ps2pdf = 'ps2pdf';

% Windows computers should have ps2pdf in their path so we don't have to
% worry about it.  Unix based computers should look for it.
missing = 1;
if ~ispc,
    [missing ps2pdf] = dos('which ps2pdf');
    if missing,
        % Here are some common paths where ps2pdf might be found.  Test them.
    	paths = {'/usr/local/bin/ps2pdf', '/usr/texbin/ps2pdf', ...
    	         '/opt/local/bin/ps2pdf', '/usr/bin/ps2pdf'};
        for location = paths,
            [missing ps2pdf] = dos(['which ' char(location)]);
            if ~missing, break; end
        end
%        if missing,
%            warning('Unable to find ps2pdf, will generate EPS files instead.');
%        end
    end
end
try
    ps2pdf = strtrim(ps2pdf);
    dos( [ps2pdf ' -dEPSCrop ' epsName ' ' pdfName] );
    delete(epsName);
catch
    warning('Unable to find ps2pdf, will generate EPS file instead.');
end

if use_file, close(fig); end


    function arrange_verticle()
        
        % Make sure the units are all 'normalized'
        plots = findobj(fig,'Type','axes')';
        for splot = plots,
            set( splot, 'Units', 'normalized' );
        end

        % The algorithm here is strange.  First we recorded a list of subplots and
        % their associated positions in the above loop.  Then, we calculate a scaling
        % factor 'scalar'.  Finally, we go through and rescale each plot.  The problem
        % occurs when scaling - we must scale both height and width.  In addition,
        % subsequent plots must also be scaled + positioned.

        max_west_legend_width = 0;
        max_east_legend_width = 0;
        
        plots = findobj(fig,'Type','axes','Tag','legend')';
        
        % Find the maximum legend width that is set to be outside the main plot axes.
        for splot = plots,
            leg_location = get( splot, 'Location' );
            if ~isempty( strfind( leg_location, 'Outside' ) ),
                legend_outer = get( splot, 'outerposition' );
                if legend_outer(1) < 0.5,
                    if max_west_legend_width < legend_outer(3),
                        max_west_legend_width = legend_outer(3);
                    end
                else
                    if max_east_legend_width < legend_outer(3),
                        max_east_legend_width = legend_outer(3);
                    end
                end
            end
        end
        
        plots = findobj(fig,'Type','axes','-not','Tag','legend')';
        for i = 1:length(plots),
            [junk, y_position] = center_axes(plots(1,i));
            plots(2,i) = y_position;
        end
        [junk, order] = sort( plots( 2,: ) );
        temp2 = [];
        for index = order,
            temp2(end+1) = plots(1,index);
        end
        plots = temp2;        
        
        % Now we have to find the largest 'TightInset' values
        inset = [];
        for splot = plots,
            inset(end+1,:) = get( splot, 'TightInset' );
        end
        if length( inset(:,1) ) > 1,
            inset = max( inset );
        end
        
        % Set the outside plot positions to use the entire area.
        startPos = 0;
        for splot = plots,
            position = [];            
            position(1) = 0;
            position(2) = startPos;  startPos = startPos + 1/length(plots);
            position(3) = 1;
            position(4) = 1/length(plots);
            
            set( splot, 'outerposition', position );
          
            position = get( splot, 'position' );            
            position(1) = max_west_legend_width + inset(1);
            position(3) = 1 - max_east_legend_width - inset(1) - inset(3);
            set( splot, 'position', position );
            
            set( splot, 'ActivePositionProperty', 'Position' );
        end
    end



    function arrange_horizontal()

        % Make sure the units are all 'normalized'
        plots = findobj(fig,'Type','axes')';
        for splot = plots,
            set( splot, 'Units', 'normalized' );
        end

        
%        plots = findobj(fig,'Type','axes','Tag','legend')';
        
        % Find the maximum legend width that is set to be outside the main plot axes.

        % Does not apply correctly to horizontally arranged plots
%         for splot = plots,
%             leg_location = get( splot, 'Location' );
%             if ~isempty( strfind( leg_location, 'Outside' ) ),
%                 legend_outer = get( splot, 'outerposition' );
%                 if legend_outer(1) < 0.5,
%                     if max_west_legend_width < legend_outer(3),
%                         max_west_legend_width = legend_outer(3);
%                     end
%                 else
%                     if max_east_legend_width < legend_outer(3),
%                         max_east_legend_width = legend_outer(3);
%                     end
%                 end
%             end
%         end
        
        plots = findobj(fig,'Type','axes','-not','Tag','legend')';
        for i = 1:length(plots),
            x_position = center_axes(plots(1,i));
            plots(2,i) = x_position;
        end
        [junk, order] = sort( plots( 2,: ) );
        temp2 = [];
        for index = order,
            temp2(end+1) = plots(1,index);
        end
        plots = temp2;        
        
        % Now we have to find the largest 'TightInset' values
        inset = [];
        for splot = plots,
            inset(end+1,:) = get( splot, 'TightInset' );
        end
        inset = max( inset );
        
        % Set the outside plot positions to use the entire area.
        startPos = 0;
        for splot = plots,
            position = [];            
            position(1) = startPos;  startPos = startPos + 1/length(plots);
            position(2) = 0;
            position(3) = 1/length(plots);
            position(4) = 1;
            set( splot, 'outerposition', position );          
            set( splot, 'ActivePositionProperty', 'Position' );
        end        
    end


    function position_legends()
        legends = findobj(fig,'Type','axes','Tag','legend')';
        l_plot = findobj(fig,'Type','axes','-not','Tag','legend')';
        
        if length(legends) ~= length(l_plot), return, end
        
        % We know the legend and associated axis will be close together.
        % Note that each axes must already have a legend
        for splot = l_plot,
            if ~strcmp(legendPosition, '')
                legend( splot, 'location', legendPosition );
            end
        end
    end

    function [centerX, centerY] = center_axes( fig_axes )
        center_position = get( fig_axes, 'position' );
        centerX = center_position(1) + center_position(3)/2;
        centerY = center_position(2) + center_position(4)/2;
    end

    function new_figure = CloneFig(inFig)
        % this program copies a figure to another figure
        % example: CloneFig(1,4) would copy Fig. 1 to Fig. 4
        % Matt Fetterman, 2009
        % pretty much taken from Matlab Technical solutions:
        % http://www.mathworks.com/support/solutions/en/data/1-1UTBOL/?solution=1-1UTBOL
%         hf1=figure(inFigNum);
        new_figure=figure();
        clf;
        compCopy(inFig,new_figure);
        
        function compCopy(op, np)
            %COMPCOPY copies a figure object represented by "op" and its % descendants to another figure "np" preserving the same hierarchy.
            
            ch = get(op, 'children');
            
            if isempty(ch), return, end

            nh = copyobj(ch,np);
            for k = 1:length(ch)
                compCopy(ch(k),nh(k));
            end
        end
    end

end
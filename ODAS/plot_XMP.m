%% plot_XMP
% Simple real-time plots of data collected with XMP instrument.
%%
% <latex>\index{Type A!plot\_XMP}</latex>
%
%%% Syntax
%   plot_XMP( fileName )
%
% * [fileName] name of the data file to display
%
%%% Description
%
% Provides a simple but effective preview of the data from a vertical 
% profile.  Plotted on a time vs count graph, the data forms descending 
% traces where each trace represents a channel.  Visualizing the data in 
% this manner lets one see what is happening to the instrument in either 
% real-time or during a previously recorded acquisition.
%
% When run, the user is asked for the figure duration.  The duration is the 
% length of the vertical axis in seconds.  The function will plot the 
% requested channels over this duration as a single figure.  When the end 
% of the figure is reached, the plot will pause before clearing the graph 
% and continuing from where it left off at the top of the graph.
%
% Data is plotted in units of counts - essentially raw values from the 
% analog to digital converter.  When plotted, the resulting traces tend to 
% overlap and are difficult to see.  To solve this problem one should 
% apply scalar and offset values to position traces to ensure they can be 
% viewed.
%
% The channels to display, along with their respective scalar and offset 
% modifiers, are controlled by variables defined near the top of this 
% function.  The section looks similar to what is shown below:
%
%    % Variables to be plotted, and their channel numbers
%    fast_vars     = {'Ax','Ay','Az','T1_dT1','T2_dT2','Sh1','Sh2','C1_dC1',
%                     'C2_dC2','Fluo','BS'};
%    fast_var_nums = [ 1  2  3  5  7  8  9  12  13  14  15 ];
%    slow_vars     = {'P','P_dP', 'Mx', 'My'};
%    slow_var_nums = [10   11      34    33];
%  
%    set(0,'Defaultaxesfontsize',10,'Defaulttextfontsize',10);
%    Ax_scale        = 1 ;    Ax_offset         =  -30000;         
%    Ay_scale        = 1 ;    Ay_offset         =  -25000;
%                Continued in function....
%

% *Version History:*
%
% * 2004-06-06 RGL
% * 2004-09-13 IG modifications for speed and appearance
% * 2005-01-18 IG? minor modifications to the size & position of the figure & axes
% * 2005-02-15 IG? no longer asks the user whether or not to continue
%   plotting at the end of each window.  User can interrupt plotting in two
%   ways: Ctrl-C or simply close the window (the routine quits if it can't
%   find the window it created). If no extension is included on the file
%   name, now assume ".p" ending.
% * 2005-02-18 IG? Change routine so that channel names and numbers are
%   provided in the parameters, and locations of each channel within the
%   matrices read in from the file are determined automatically.
% * 2005-02-21 IG? Modified approach to plotting, plotting time for a 100s
%   window cut roughly in half. May be slightly more memory-intensitve,
%   however, and assumes that OpenGL graphics acceleration is available.
% * 2005-07-12 IG? Now assume two microconductivity channels as the "default"
%   setup for an instrument.
% * 2005-07-15 IG? Added clause to deal with variables in the default setup
%   that are not actually present in the data file
% * 2006-05-14 RGL  Added little endian flag to fileopen. Corrected pause at
%   end of file while reading data.
% * 2009-03-06 RGL changed fopen to fopen_odas.
% * 2010-01-15 AWS support odas V6 and up
% * 2010-10-05 AWS added doxygen tags
% * 2012-01-16 RGL special version made for Kats at JAMSTEC
% * 2012-04-11 WID changed inifile_with_instring calls to setupstr
% * 2012-09-09 WID updated documentation
% * 2013-02-26 WID merge of changes from Rolf.  Allow matrix / rate of any size.
%                  Previously would only work with expected values.


function plot_XMP(fileName)

% parsing of command line parameters
if nargin < 1; fileName    = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get Input from the User
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(fileName); fileName = getDataFileName; end;
if (fileName == 'q'); return;end

max_plot_length_in_records = input('How many seconds of data per plot? (default=100s) ');
if isempty(max_plot_length_in_records), max_plot_length_in_records=100; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization of constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables to be plotted, and their channel numbers
fast_vars     = {'Ax','Ay','T1_dT1','Sh1','Sh2'}; 
fast_var_nums = [ 1    2      5       8     9    ];
slow_vars     = {'P','P_dP', 'Pitch' 'V_Bat'};
slow_var_nums = [10   11       13      32   ];

% Plotting parameters. Offsets etc. may be listed for as many potential
% plotting variables as are desired; routine will only plot those listed
% above.
% Customized for XMPs used by Kats at JAMSTEC
set(0,'Defaultaxesfontsize',10,'Defaulttextfontsize',10);
Ax_scale        = 1 ;    Ax_offset         =  -30000;         
Ay_scale        = 1 ;    Ay_offset         =  -25000;
T1_dT1_scale    = 1 ;    T1_dT1_offset     =  -10000;
Sh1_scale       = 1 ;    Sh1_offset        =   0;
Sh2_scale       = 1 ;    Sh2_offset        =   10000;
P_dP_scale      = 4 ;    P_dP_offset       =  -50000;
P_scale         = 4 ;    P_offset          =  -50000;
Pitch_scale     = 1;     Pitch_offset      =   30000;
V_Bat_scale     = 1;     V_Bat_offset      =   20000;

% Factors by which to decimate fast & slow data for plotting (decimation is
% used to speed up plotting of long records)
if max_plot_length_in_records>=500
    decimate_fast = 256;
    decimate_slow = 32;
elseif max_plot_length_in_records>=250
    decimate_fast = 128;
    decimate_slow = 16;
elseif max_plot_length_in_records>=100
	decimate_fast = 64;         
	decimate_slow = 8;
elseif max_plot_length_in_records>=50
    decimate_fast = 8;
    decimate_slow=1;
else
    decimate_fast=1;
    decimate_slow=1;
end
figurePos = [0.005 0.03 0.65 0.89];  % Figure position
figureUnits = 'normalized';
axesPos = [0.07 0.06 0.8 0.85];
xlims = [-60000 40000];
if length(fast_vars)<=8, cmap = [0 0 0.5; 0 0 1; 0 1 1; 0 0.5 0; 0 1 0; 0.5 0 0; 1 0 0; 1 0 1; 0 0 0; 0.5 0.5 0.5; 1 1 0];
else cmap = [0 0 0.5; 0 0 1; 0 1 1; 0 0.5 0; 0 1 0; 0.5 0 0; 1 0 0; 1 0 1; 1 1 0; 0 0 0; 0.5 0.5 0.5]; end

% ODAS parameters
size_of_integer =  2;% size of integers from data files in bytes
header_size_i   = 18;% index of header_size in header record
block_size_i    = 19;% index of block_size in header record

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Open the file, get required information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fileOpen(fileName); % open file to get some info out of its header

% extract some record parameters  
% assume that header has at least 64 2-byte words
HD = fread(fid, 64, 'ushort');	% read header
frewind(fid);			% back to the beginning of the file
header_size_in_bytes = HD(header_size_i); % This is the actual size of the header
header_size = header_size_in_bytes / size_of_integer;
data_record_size_in_bytes = HD(block_size_i)-header_size_in_bytes;
data_record_size = data_record_size_in_bytes / size_of_integer;
outbound_clock = HD(21) + HD(22)/1000; % Outbound address clock in Hz
header_version = bitshift(HD(11), -8) + bitand(HD(11), 255) /1000;
if header_version >= 6
    first_record_size_in_bytes = header_size_in_bytes + HD(12);
else
    first_record_size_in_bytes = data_record_size_in_bytes + header_size_in_bytes;
end

% ODAS4IR has a bug. The clock frequency is not written into the header
% until the second header.

fseek (fid, first_record_size_in_bytes,'bof');
second_header = fread(fid, 64, 'ushort');
frewind(fid); % back to beginning of file
outbound_clock = second_header(21) + second_header(22)/1000; % Outbound address clock in Hz
record_duration = data_record_size / outbound_clock; % the length of a record in seconds

date_string =[num2str(HD(4)) '-' num2str(HD(5)) '-' num2str(HD(6)) '   ' num2str(HD(7)) ':' ...
        num2str(HD(8)) ':' num2str(HD(9)) '.' num2str(HD(10))];
% Need to get sampling rate

% load channel info by using the address matrix in record number zero of the data file
% compute sizes of setup matrix, total number of channels, etc.
[rows, cols, no_slow_cols, slow_ch, fast_ch] = load_ch_setup(fid);
no_of_ch = length(slow_ch) + length(fast_ch); % Total number of data channels in the file	
fast_points_per_record = data_record_size / cols; % number of data points of a fast channel in a single record.
slow_points_per_record = fast_points_per_record / rows; % Number of points of a slow channel in a single record

max_length_of_fast_channels = max_plot_length_in_records * fast_points_per_record;% 
max_length_of_slow_channels = max_plot_length_in_records * slow_points_per_record;% 
for ii=1:length(fast_vars)      % Figure out positions of variables in fast & slow matrices
    eval([fast_vars{ii} '_ch = find(fast_ch==' num2str(fast_var_nums(ii)) ');']);
end
for ii=1:length(slow_vars)
    eval([slow_vars{ii} '_ch = find(slow_ch==' num2str(slow_var_nums(ii)) ');']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading & Plotting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hF=figure(1);
set(hF,'units',figureUnits,'position',figurePos,'userdata','XMP_fig','renderer','openGL');
clf
h =axes('position',axesPos);

for ii=1:length(fast_vars)          % Prepare the text for the legend
    junk1 = eval([fast_vars{ii} '_scale']); 
    if junk1==1, junk1=''; else junk1=[num2str(junk1) '*']; end
    junk2 = eval([fast_vars{ii} '_offset']); 
    if junk2<0, junk2=num2str(junk2); elseif junk2==0, junk2=''; else junk2=['+' num2str(junk2)]; end
    leg_text{ii} = [junk1 strrep(fast_vars{ii},'_','\_') junk2];
end
for ii=1:length(slow_vars)
    junk1 = eval([slow_vars{ii} '_scale']); 
    if junk1==1, junk1=''; else junk1=[num2str(junk1) '*']; end
    junk2 = eval([slow_vars{ii} '_offset']); 
    if junk2<0, junk2=num2str(junk2); elseif junk2==0, junk2=''; else junk2=['+' num2str(junk2)]; end
    leg_text{ii+length(fast_vars)} = [junk1 strrep(slow_vars{ii},'_','\_') junk2];
end
t_f = (0:max_length_of_fast_channels -1)'/ fast_points_per_record;% time vector for plotting of fast channels
t_s = (0:max_length_of_slow_channels -1)'/ slow_points_per_record;% time vector for plotting of slow channels

fast_data = zeros(fast_points_per_record,length(fast_vars))*NaN;   % pre-assign matrices
slow_data = zeros(slow_points_per_record,length(slow_vars))*NaN;

status=fseek(fid,first_record_size_in_bytes,'bof'); % move to the beginning of the first real data record
if status==-1
    close(1)
    error(['Error trying to move to the beginning of the first real data record: ' ferror(fid)]);
end

d_buffer        = ones (data_record_size,1)*NaN;
h_buffer        = ones (header_size,1)*NaN;

we_are_not_done = 1;    % flag for stopping the reading and plotting
fast_index = 1;
slow_index = 1;
record_counter = 0;
select_slow = 0;
select_fast = 0;

figure(1); cla
ylims = [0 max_plot_length_in_records];
set(h,'clipping', 'off','ylim',ylims,'ydir','rev','xlim',xlims);
set(h,'ColorOrder',cmap);
ylabel('\it t \rm [s]','fontsize',11); 
title([strrep(fileName,'_','\_') ';  ' date_string ' UT'],'fontsize',12,'fontweight','bold'); hold on; grid on
new_plot=1;
Y_all = zeros(floor(max_length_of_fast_channels/decimate_fast),length(fast_vars))*NaN;  % Initialize plotting matrices
t_all = zeros(floor(max_length_of_fast_channels/decimate_fast),1)*NaN;
Y2_all = zeros(floor(max_length_of_slow_channels/decimate_slow),length(slow_vars))*NaN;
t2_all = zeros(floor(max_length_of_slow_channels/decimate_slow),1)*NaN;
while (we_are_not_done)
    [h_buffer, count] = fread (fid, header_size,'short');
    if count ~= header_size, % we have a read failure.
        fseek(fid,-count,'cof');% Move back to where we were in the file.
        pause(2*record_duration);% Pause for the duration of 2 records and try again. 

        [h_buffer, count] = fread (fid, header_size,'short');
        if count ~= header_size, % we have the second read failure.
            fseek(fid,-count,'cof');% Move back to where we were in the file.
            pause(2*record_duration);% Pause for the duration of 2 records and try again.
            
                [h_buffer, count] = fread (fid, header_size,'short');
                if count ~= header_size, % we have the third read failure. 
                    % This must be the hard end of file. Give up on this file.
                    we_are_not_done = 0; 
                    'Pause for header exceeded the durtion of 4 records'%
                end
        end
    end

    if (we_are_not_done)
        [d_buffer, count] = fread (fid, data_record_size,'short'); %Read a record of data
        if count ~= data_record_size;
            fseek(fid,-count,'cof');% Move back to where we were in the file.
            pause(2*record_duration);% Pause for 2 record durations and try again.

            [d_buffer, count] = fread (fid, data_record_size,'short'); 
            if count ~= data_record_size, 
                 % Second read attempt failed. We have a serious error because
                 % the header was read but there is no data after it.                 
                we_are_not_done = 0;
                'No data behind the header even after pausing for duration of 2 records' 
            end
        end
        record_counter = record_counter + 1;
    end

    if (we_are_not_done)
        [fast_data_all, slow_data_all] = demultiplex (d_buffer, h_buffer); % Put data into matrices
        for ii=1:length(fast_vars)
            junk = eval(['fast_data_all(:,' fast_vars{ii} '_ch)*' fast_vars{ii} '_scale +' fast_vars{ii} '_offset;']);
            if ~isempty(junk)
                fast_data(:,ii) = junk;
            else
                fast_data(:,ii) = ones(size(fast_data_all,1),1)*NaN;
            end
        end
       for ii=1:length(slow_vars)
            junk = eval(['slow_data_all(:,' slow_vars{ii} '_ch)*' slow_vars{ii} '_scale +' slow_vars{ii} '_offset;']);
            if ~isempty(junk)
                junk = junk(:,1); % in case of multiple samples in slow channels
                slow_data(:,ii) = junk;
            else
                slow_data(:,ii) = ones(size(slow_data_all,1),1)*NaN;
            end
        end
        clear fast_data_all slow_data_all
        ii=fast_index:fast_index+fast_points_per_record - 1;
        ii2=slow_index:slow_index+slow_points_per_record - 1;
        fast_index = fast_index + fast_points_per_record;
        slow_index = slow_index + slow_points_per_record;
        Y = decimate_me(fast_data,decimate_fast); % reduce the number of points
        t = decimate_me (t_f(ii),decimate_fast); % also for time vector

        select_fast = select_fast(end) + (1:size(Y,1))';
        
        Y_all(select_fast,:)=Y; t_all(select_fast)=t;             % Add data to the plotting matrices
        
        Y2 = decimate_me(slow_data,decimate_slow);
        t2=decimate_me(t_s(ii2),decimate_slow);
        select_slow = select_slow(end) + (1:size(Y2,1))';
        
%        ii = (record_counter-1)*fast_points_per_record/decimate_fast+1:record_counter*fast_points_per_record/decimate_fast

        
%        ii = (record_counter-1)*slow_points_per_record/decimate_slow+1:record_counter*slow_points_per_record/decimate_slow

        Y2_all(select_slow,:)=Y2; t2_all(select_slow)=t2;

        % setup timer to limit our graphing commands to 1 per second.
        if ~exist('lasttime', 'var'),
            lasttime = clock - 2;         % Set to trigger on the first loop
            draw_refresh_count = 0;
        end
        
%        if ~isempty(findobj('userdata','XMP_fig')), % Update the plot
            if record_counter==1, hh1=plot(Y_all,t_all); hh2=plot(Y2_all,t2_all);
            else
              if draw_refresh_count == 0 && etime(clock, lasttime) > 1,
                for ii=1:size(Y_all,2), set(hh1(ii),'xdata',Y_all(:,ii),'ydata',t_all); end
                for ii=1:size(Y2_all,2), set(hh2(ii),'xdata',Y2_all(:,ii),'ydata',t2_all); end
              end
            end
%        else we_are_not_done=0; return; end
        
        if new_plot         % Only re-create the legend when starting a new plot
            legend(leg_text,-1);
            new_plot=0;
        end
        if draw_refresh_count == 0,
            draw_refresh_count = 0;
            if etime(clock, lasttime) > 1,
                lasttime = clock;
                drawnow
            end
        else
            draw_refresh_count = draw_refresh_count - 1;
        end
    else
        for ii=1:size(Y_all,2), set(hh1(ii),'xdata',Y_all(:,ii),'ydata',t_all); end
        for ii=1:size(Y2_all,2), set(hh2(ii),'xdata',Y2_all(:,ii),'ydata',t2_all); end
        drawnow
    end
    
    if record_counter == max_plot_length_in_records         % When done a plot, clean up the figure and prepare to plot some more
        for ii=1:size(Y_all,2), set(hh1(ii),'xdata',Y_all(:,ii),'ydata',t_all); end
        for ii=1:size(Y2_all,2), set(hh2(ii),'xdata',Y2_all(:,ii),'ydata',t2_all); end
        drawnow
        pause(1)
        figure(1); cla
        Y_all(:,:)=NaN; t_all(:,:)=NaN; Y2_all(:,:)=NaN; t2_all(:,:)=NaN;
        ylims = ylims+max_plot_length_in_records;
        set(h,'clipping', 'off','ylim',ylims,'ydir','rev','xlim',xlims);
    % 	cmap = flipud(ColorCube(length(slow_vars)+length(fast_vars)+2)); cmap=cmap(2:end,:);
        set(h,'ColorOrder',cmap);
        ylabel('\it t \rm [s]','fontsize',11); 
        title([strrep(fileName,'_','\_') ';  ' date_string],'fontsize',12,'fontweight','bold'); hold on; grid on
        t_f = t_f + max_length_of_fast_channels/fast_points_per_record;
        t_s = t_s + max_length_of_slow_channels/slow_points_per_record;
        fast_index = 1;
        slow_index = 1;
        record_counter = 0;
        select_slow = 0;
        select_fast = 0;
    end    
end

%******************************************************************
function Y = decimate_me(X,R)
%
% function to speed up plotting by showing only max and min
% of successive data ensembles R long
% R is the reduction ratio, a positive integer and
% The number of rows in Y is R times smaller than the rows in X.
% X and Y have the same number of columns 
% 
% RGL 2013-02-07

R = floor(abs(R));% In case R is not a positive integer

% decimate vectors by factor R
len = size(X,1);

X = X(1:2*R*floor(len/(2*R)),:); % truncate, if necessary
len = size(X,1); % The (possibly) truncated number of rows of X

junk = reshape(X, 2*R, len/(2*R), size(X,2));
Y = reshape([min(junk); max(junk)], len/R, size(X,2)); 

 
%******************************************************
function [fast_data, slow_data]=demultiplex(data,header)
% [fast_data, slow_data]=demultiplex(data,header)
%
% function to demultiplex data collected with the ODAS system.
% data is the vector of multiplexed data as recorded by ODAS.
% The header is used to figure out the number of slow and fast columns and the number of rows in 
% basic matrix used to multiplex the data.

if (~(nargin ==2 )); error ('input must have 2 arguments'),end;

fast_columns = header(29);% Matrix information is located in the header
slow_columns = header(30);
rows =         header(31);

columns = fast_columns + slow_columns ; % define total number of columns in data
data = reshape(data,columns,length(data)/columns);
data = data';
fast_data = data(:,slow_columns+1:columns);% extract the fast channels

n_slow = rows*slow_columns; % Number of slow channels
if (slow_columns ~= 0)
        slow_data = data(:,1:slow_columns); % extract slow columns
        slow_data = slow_data';
        slow_data = slow_data(:);
        slow_data = reshape(slow_data,n_slow,length(slow_data)/n_slow);
        slow_data = slow_data'; % put slow vectors into columns
else
        slow_data = [];
end

%******************************************************************************
function fileName = getDataFileName;
% fileName = getDataFileName returns the file name of the data file.
% Preforms error checking.
%
% Fab, March 1998.
% Modified by RGL to offer the latest '*.p' file for opening.
% 2004-06--6


fileName='';
while isempty(fileName); 
    test_string = get_latest_file;
   fileName = input(['Enter data file name (default: ' test_string ' , ''q'' to quit): '], 's');
   if strcmp(fileName, 'q')
      fclose('all');return;
   elseif isempty(fileName)
      fileName = test_string;
   end
   if ~exist(fileName, 'file');
      fprintf('Can''t open file %s! Try again.\n', fileName);
      fileName = '';
   end
end


%******************************************************************************
function fid = fileOpen(fileName);
% fid = fileOpen(fileName) returns the file ID for the file fileName.
% Preforms error checking.
%
% Fab, March 1998.

fid = -1;
[fid, error_message] = fopen_odas(fileName, 'r');
if ~isempty(error_message), error_message, end;
if fid == -1,
   error(sprintf('Error opening file %s !\n', fileName));
   return
end

%******************************************************************************
function [nRow, nCol, nSlowCol, slowCh, fastCh] = load_ch_setup(fid)
% [nRow, nCol, nSlowCol, slowCh, fastCh] = load_ch_setup(fid) 
% loads the channel matrix from the file fileName.
% Returns the number of rows in the matrix, the number of columns in the matrix,
% the number of slow columns, a vector containing the slow channels,
% and a vector containing the fast channels.
% 
% Replaces load_ch_setup.m written by L. Zhang.
% Fab, March 1998.
% AWS - 2010-01-14 changes for odas v6 and up

fseek(fid,0,'bof'); % Rewind to the beginning of the file
header = fread(fid, 64, 'ushort');%read the header
header_version = bitshift(header(11), -8) + bitand(header(11), 255) /1000;

if header_version >= 6
    nFast = header(29);
    nSlowCol = header(30);
    nRow  = header(31);
    nCol = nSlowCol + nFast;
    matrixSize = nRow * nCol; %use this to check what we get out of the setup file string
    setupfilestr = char(fread(fid, header(12), 'char'));
    
    if isempty(setupfilestr)
        error('failed to extract setup file string from first data record');
    end
    setupfilestr = setupfilestr';
    
    cfg = setupstr(setupfilestr);
    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    matrix = [];
    for row = rows
      values = textscan(row{1}, '%d16');
      matrix = vertcat(matrix, values{1}');
    end
    
    slowCh = matrix(:, 1:nSlowCol)';
    slowCh = slowCh(:);
    if length(slowCh) ~= nSlowCol*nRow
        error('Error building channel matrix: number of slow channels does not agree.');
    end

    fastCh = matrix(1, nSlowCol+1:nCol)';
    if length(fastCh) ~= nFast
        error('Error building channel matrix: number of fast channels does not aggree.');
    end

else

    nFast = header(29);% extract info
    nSlowCol = header(30);
    nRow = header(31);

    % rebuild the channel matrix
    nCol = nSlowCol + nFast;
    matrixSize = nRow * (nSlowCol+nFast);
    matrix = fread(fid, matrixSize, 'ushort');
    matrix = reshape(matrix, nCol, matrixSize/nCol)';


    % build the slowCh vector
    slowCh = matrix(:, 1:nSlowCol)';
    slowCh = slowCh(:);
    if length(slowCh) ~= nSlowCol*nRow,
        error('Error building channel matrix: number of slow channels does not agree.');
    end

    %build the fastCh vector
    fastCh = matrix(1, nSlowCol+1:nCol)';
    if length(fastCh) ~= nFast,
       error('Error building channel matrix: number of fast channels does not aggree.');
    end
end

return

%********************************************************************
function     test_string = get_latest_file();
%
% function to get the name of the "*.p" file in the current directory with
% the latest date. The assumption is that the user wants to plot the
% latest data file.
% RGL 2004-06-06

D = dir('*.p');
if (isempty(D)); test_string=[]; return; end% There are no *.p files in this directory
test_string = D(1).name; % use the first file in the list of *.p files
creation_date = datenum(D(1).date);% The date and time of its creation
if size(D,1) > 1 ; % more than 1 file was found, so test which is the latest
    for k = 2:size(D,1)
        if creation_date < datenum(D(k).date)
            creation_date = datenum(D(k).date); % This file is newer so use it.
            test_string = D(k).name;
        end
    end
end

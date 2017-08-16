%% read_odas
% Convert an ODAS binary file to .mat format
%%
% <latex>\index{Type A!read\_odas}</latex>
%
%%% Syntax
%   [var_list, outname] = read_odas( fname, ch_nums, ch_names )
%
% * [fname] Name of ODAS binary file. Optional, user prompted if omitted.
% * [ch_nums] Vector of channel numbers to convert. Optional, all channels
%       selected by default.
% * [ch_names] Vector of names to assign the channels provided by 'ch_num'.
%       Optional, default names used if not provided.
% * []
% * [var_list] List of variables within the resulting .mat file.
% * [outname] Name of the resulting .mat file.
%
%%% Description
%
% Convert an ODAS binary file (.p) into a MATLAB file (.mat).  All requested 
% channels are extracted from the input file and saved within a .mat file.  The
% resulting .mat file can be loaded directly into Matlab.
%
% If ch_nums is ommited, all channels found within the input file are 
% extracted. Specific channels can be extracted by including ch_nums - a vector
% of numbers matching the requested channel numbers.
%
% The extracted channels are given default names. Perferred names can be 
% provided by including a ch_names argument. The ch_names argument should be a 
% vector containing names to assign to those channels identified by ch_nums.
% If ch_names is included, ch_names and ch_nums should be vectors of the same
% length. ODAS v6 data files do not require this argument as the channel names 
% are embedded inside the file.
% 
% If default channel names are used and the function can not guess what the 
% correct name should be, a standard name is generated consisting of 'ch' 
% and the channel number.  For example, channel 239 would be named 'ch239'.
%
% The default .mat file name is the input file name, fname, with the extension
% changed to '.mat'.  If this file already exists, the user is queried for an
% alternate file name.  When the output argument outname is used, a temporary 
% name that differs from the default name is used for the .mat file. This is 
% useful when functions need to access the raw data and will delete the 
% generated .mat file when finished.
%
%%% Examples
%
%    >> read_odas;
%
% Extract all data from a .p file and store inside a newly generated .mat file.
% The input file name will be requested from the user.
%
%    >> my_vars = read_odas('my_file.p');
%
% Extract all data from 'my_file.p' and store inside 'my_file.mat'.  Returns 
% a list of the extracted channels inside 'my_vars'.
%
%    >> my_vars = read_odas('my_file.p', [1:4], {'gnd1','dw1','SBT1E','SBT1O'});
%
% Function used on a legacy data file (pre-v6 header).  Extract channels 
% 1 to 4 from 'my_file.p' and store in 'my_file.mat'.  The extracted channels
% are named 'gnd1', 'dw1', 'SBT1E', 'SBT1O' and are returned inside 'my_vars'.
%
%%% Legacy ODAS Notes (pre v6)
%
% # Default channel numbers and names are defined at the start of this function 
% and can be modified for your instrument.
% # Split channels are used when a sensor generates data as 32bit words. These 
% words are addressed using two sequential 16bit channels. When read, these 
% channels should be addressed individually with 'E' and 'O' suffixes appended 
% to their names.  When found, read_odas will internally combine the channels 
% and return a single value. This behaviour can be turned off by modifying the 
% 'combine_channels' flag inside the function.

% *Version History:*
%
% * 2005-02-01 (IG) Based in part on misc. OTL/RGL routines, including
%   plot_tomi, demultiplex, read_tomi, convert_header, convert_mat_slow,
%   convert_mat_fast, and others
% * 2005-05-02 (IG) added second C_dC channel to default channel number/name
%   setup to deal with Peter Winsor's instrument
% * 2006-04-17 (RGL) added 1e-12 to data so that all variables are saved
%   as true 8-byte floats. Also forced a -v6 save to mat-file
% * 2007-11-05 (RGL) added SBE41F oxygen sensor to channels 48 and 49.
% * 2009-03-06 (RGL) changed fopen to fopen_odas.
% * 2009-12-22 (RGL) changed to demultiplexing record-by-record for better
%   memory efficiency.
% * 2010-01-15 (AWS) added support for odas v6 and up
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% * 2012-07-19 (WID) speed improvements - reduced load from waitbar
% * 2012-10-23 (WID) documentation update and small logic change when opening 
%                    the .mat file.
% * 2012-11-23 (WID) octave compatability change when setting output file name

function [var_list, outname] = read_odas(fname,ch_nums,ch_names)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% File handling flags and channel setup                   %%%%%
%%%%%% This section may be changed to reflect instrument setup %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default channel setup: ch_nums is a vector of channel numbers that may be
% used; ch_names is a string array of corresponding channel names (these
% are the names that the variables will bear).
if nargin<3
    
% The following two sets of ch_nums and ch_names are only relevant for
% for legacy, v1 header data files. Data files with v6 headers will
% automatically adapt.

% Default setup for most instruments
    ch_nums = [0:15 16:19 32:37 48:49 255];
    ch_names = {'gnd', 'Ax', 'Ay', 'Az', 'T1', 'T1_dT1', 'T2', 'T2_dT2', ...
        'sh1', 'sh2', 'P', 'P_dP', 'C1_dC1', 'C2_dC2', 'ch14', 'ch15', ...
        'sbtE', 'sbtO', 'sbcE', 'sbcO', ...
        'Mz', 'My', 'Mx', 'Ux', 'Uy', 'Vbat', ...
        'O2E', 'O2O', ...
        'sp_char'};

%{  
% Default setup for VMP-200 instruments
    ch_nums = [0:15 16:19 32:33 40:42 255];
    ch_names = {'gnd', 'Ax', 'Ay', 'Az', 'T1', 'T1_dT1', 'T2', 'T2_dT2', ...
        'sh1', 'sh2', 'P', 'P_dP', 'PV', 'ch13', 'ch14', 'ch15', ...
        'sbtE', 'sbtO', 'sbcE', 'sbcO', ...
        'V_Bat', 'C1_dC1', ...
        'Incl_X', 'Incl_Y', 'Incl_T', ...
        'sp_char'};
%}    
end

% File handling flags
combine_channels = 1;       % set to 1 to combine odd/even channels, 0 to leave as is

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Check file name, check for a .mat file, and open the .p file %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask the user for a file name if not provided as an input variable
if nargin<1
    fname = input(['Enter data file name (default: ' get_latest_file '): '],'s');
    if isempty(fname), fname = get_latest_file; end
end

[P,N,E] = file_with_ext( fname, {'','.p','.P'} );
if isempty(N), error('Unable to find file %s\n',fname); end

outname = N;
if nargout == 2
    outname = [N '_tmp.mat'];
elseif exist([outname '.mat'], 'file')
	sub_name = input(['.mat file already exists. \nEnter a new name to ' ...
                      'save to a different file, "Enter" to overwrite, ' ...
                      'or "q" to quit: '], 's');

    if strcmpi(sub_name, 'q'), return; end
    if ~isempty(sub_name), outname = sub_name; end
end

% Ensure the output name ends in '.mat'.  Required for octave.
[jP, jN, jE] = fileparts(outname);
outname = [jN '.mat'];

% Open the file
[fid, error_message] = fopen_odas([N,E],'r');
if ~isempty(error_message), disp(error_message); end;
if fid<3
   error('Unable to open file %s\n',[N,E]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Constants & parameters                                   %%%%%%
%%%%%% These will not generally need to be changed by the user. %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
bytes_per_word      =  2; % Assume 2-byte integers
header_version_i    = 11; % this contains header version, always 1 prior to odas v6
setupfile_size_i    = 12; % odas v6 and higher only: size of ini string in first record
header_size_i       = 18; % Index of header_size in header record
block_size_i        = 19; % Index of block_size in header record
header_size_0       = 64; % Predicted size of the header in each block (integers)
fast_cols_i         = 29; % Header indices for the number of fast & slow columns
slow_cols_i         = 30;
rows_i              = 31; % Header index for the number of rows
clock_whole_i       = 21; %Index to integer part of clock frequency
clock_fraction_i    = 22; %Index to fractional part of clock frequency
%a = realmin*1e100;            % Factor needed to force MATLAB to convert header data into 8-byte floats
a = 1e-10;            % Factor needed to force MATLAB to convert header data into 8-byte floats
B = 2^16;                 % Multiplier for 4-byte integers

% Parameters read in from the first header in the file
fseek(fid,0,'bof');
HD = fread(fid,header_size_0,'ushort');
header_size = HD(header_size_i)/bytes_per_word;
setupfile_size = HD(setupfile_size_i);
record_size = HD(block_size_i)/bytes_per_word;
data_size = record_size - header_size;
fast_cols = HD(fast_cols_i);
slow_cols = HD(slow_cols_i);
n_cols = fast_cols + slow_cols;
n_rows = HD(rows_i);
rows_per_record = data_size/n_cols; % rows in a record
slow_samples_per_record = rows_per_record / n_rows; % samples in a slow channel per record
n_slow = n_rows*slow_cols; % number of slow channels
f_clock = HD(clock_whole_i) + HD(clock_fraction_i)/1000;
fs_fast = f_clock / n_cols; % the sampling rate of the fast channels
fs_slow = fs_fast / n_rows; % the sampling rate of the slow channels
fseek(fid,0,'eof');
filesize = ftell(fid);

%MSB has major version, LSB has minor version
header_version = bitshift(HD(header_version_i), -8) + bitand(HD(header_version_i), 255) /1000;

if (header_version >= 6)
    %odas v6 and higher stores setup file in first record, which is variable size
    first_record_size = header_size*bytes_per_word + setupfile_size;
else
    first_record_size = record_size*bytes_per_word;
end

n_records = (filesize - first_record_size)/(record_size*bytes_per_word);

if (n_records-floor(n_records))~=0
    n_records = floor(n_records);
    warning(['The file ' P N E ' does not contain an integer number of records']);
end
if n_records <= 1
    error(['File ' P N E ' contains no data'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Read the data file %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the matrices
header = zeros(n_records,header_size)*NaN;    % Note that the first block doesn't contain any data
if fast_cols ~= 0
    fast_data = zeros(rows_per_record*n_records,fast_cols);
end
if slow_cols ~= 0
    slow_data = zeros(slow_samples_per_record*n_records,n_slow);
end

fseek(fid,header_size*bytes_per_word,'bof');

isXMP = 0;

if header_version >= 6
    clear ch_names ch_nums;
    setupfilestr = fread(fid,setupfile_size,'*char*1')';
    if (isempty(setupfilestr))
        error('failed to extract setup file string from first record');
    end
    
    % Load the matrix
    cfg = setupstr(setupfilestr);
    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    ch_matrix = [];

    for row = rows
      values = textscan(row{1}, '%d16');
      ch_matrix = vertcat(ch_matrix, values{1}');
    end

    sections = setupstr(cfg, '');
 
    if ~isempty(setupstr(cfg, '', 'xmp')), isXMP=1; else isXMP=0; end

    if (isXMP)
        index_to_serial_number = find(ch_matrix' == 254);
        bad_buffer_flag_i = 15;

        if isempty(index_to_serial_number)
            error('channel number 254 is not in address matrix');
        end

        tmphbuf =fread(fid,header_size,'ushort');
        tmpdbuf =fread(fid,data_size,'ushort');

        while tmphbuf(bad_buffer_flag_i) == 1
            display('skipping bad buffer...');
            tmphbuf =fread(fid,header_size,'ushort');
            tmpdbuf =fread(fid,data_size,'ushort');
        end

        scan_length = n_rows*n_cols;
        tmpdbuf = reshape(tmpdbuf,scan_length,data_size/scan_length);

        instrument_serial_number_section = ['xmp_' num2str(tmpdbuf(index_to_serial_number),'%05d')];
 
        % find section that contains parameters for our instrument
        if isempty(strcmp(instrument_serial_number_section,sections))
            error(['could not find the following section in the xmp configuration file: ' instrument_serial_number_section]);
        end

        ch_names = { 'Ax','Ay','P','P_dP','gnd','V_Bat','T','T_dT','Pitch', 'sh1','sh2','spare'};
        for i = 1:length(ch_names)
            try
                ch_nums(i) = str2double(setupstr(cfg, instrument_serial_number_section, [ch_names{i} '_id']));
            catch err
                error([ch_names{i} ' id is missing from configuration file']);
            end
        end

        clear tmphbuf tmpdbuf i;
      

        fseek(fid, first_record_size,'bof');

    else
        all_unique_channels = unique(sort(ch_matrix(:)));
        nch = 1;

        for section_name = sections

            %look in each section whether it contains id, id_even or id_odd. If so, then
            %we are dealing with a channel section.            
            tmpid = str2double(setupstr(cfg, section_name, 'id'));
            tmpide = str2double(setupstr(cfg, section_name, 'id_even'));
            tmpido = str2double(setupstr(cfg, section_name, 'id_odd'));
            
            if isempty(setupstr(cfg, section_name, 'id.*'))
              continue;
            end
                       
            %use the parameter 'name' from the channel section, or fall back to the section name.
            namestr = char(setupstr(cfg, section_name, 'name'));
            if isempty(namestr)
                namestr = char(section_name);
            end        
            
            for channel = all_unique_channels'
                if tmpid == channel
                    fprintf(1, '     channel: %2d = %s\n', tmpid, namestr);
                    ch_nums(nch) = tmpid;
                    ch_names{nch} = namestr;
                    nch = nch + 1;
                    break; %normal channel found, no need to look further
                end
                if tmpide == channel
                    fprintf(1, 'even channel: %2d = %s\n', tmpid, namestr);
                    ch_nums(nch) = tmpide;
                    tmpstr = strcat(namestr,'E');
                    ch_names{nch} = tmpstr;
                    nch = nch + 1;
                end
                if tmpido == channel
                    fprintf(1, ' odd channel: %2d = %s\n', tmpid, namestr);
                    ch_nums(nch) = tmpido;
                    tmpstr = strcat(namestr,'O');
                    ch_names{nch} = tmpstr;
                    nch = nch + 1;  
                end
            end
        end        
    end

else
    d_buffer = fread(fid,data_size,'short');
    ch_matrix = reshape(d_buffer(1:n_cols*n_rows),n_cols,n_rows)';
end

%the file pointer should be at the beginning of the second record

if fast_cols
    ch_fast = ch_matrix(:,slow_cols+1:end);
    ch_fast = ch_fast(1,:);
else
    ch_fast=[];
end
if slow_cols
    ch_slow = ch_matrix(:,1:slow_cols)';
    ch_slow = ch_slow(:);
else
    ch_slow=[];
end

if header_version < 6
    %hack for bug in old odas where clock is not stored in first header
    %this should also be run if v1 data file gets patched with v6 config
    %file
    HD = fread(fid,header_size,'ushort');
    f_clock = HD(clock_whole_i) + HD(clock_fraction_i)/1000;
    fs_fast = f_clock / n_cols; % the sampling rate of the fast channels
    fs_slow = fs_fast / n_rows; % the sampling rate of the slow channels
    fseek(fid, first_record_size,'bof');
end


clear ch_matrix header_size_i block_size_i header_size_0 fast_cols_i slow_cols_i rows_i HD

% Read the file
h_wb=waitbar(0,'Reading the data');
wb_count = 0;
wb_clock = clock;

for index = 1:n_records
    h_buffer = fread(fid,header_size,'ushort'); % Read in a header; if fail, trim what we failed to read
    d_buffer = fread(fid,data_size,'short');    % Read in a data record; if fail, trim what we failed to read

    header(index,:) = h_buffer'; % Add data to the growing "header" & "data" variables
    junk = reshape(d_buffer, n_cols, data_size/n_cols); junk = junk' + a; % demultiplex
            %into standard matrix and add small number to force DP floats
    fast_range = 1 + (index-1)*rows_per_record : index*rows_per_record; % range
                        % of rows to be filled for fast channels
    slow_range = 1 + (index-1)*slow_samples_per_record : index*slow_samples_per_record;
                    % range of rows to be filled for slow channels

    if fast_cols ~=0
        fast_data(fast_range,:) = junk(:,1+slow_cols:n_cols);
    end
    if slow_cols ~= 0
        junk_slow = junk(:,1:slow_cols); 
        if slow_cols > 1
            junk_slow = junk_slow'; junk_slow = junk_slow(:);
        end
        slow_data(slow_range,:) = reshape(junk_slow, n_slow, slow_samples_per_record)';
    end
    
    wb_count = wb_count + 1;
    if wb_count == 500,                      % prevent to many calls to 'clock'
        wb_count = 0;                  
        if etime(clock, wb_clock) > 1/5,     % maximum refresh rate of 5 fps
            wb_clock = clock;
            waitbar(index/(n_records),h_wb);
        end
    end
    
end

close(h_wb);
clear h_buffer d_buffer
fclose(fid);



% Assign to appropriate variable names
n_var = 1; var_list={};
for ii=1:length(ch_slow)
    jj=find(ch_nums==ch_slow(ii));
    if ~isempty(jj)
        junk = char(ch_names{jj});
    else
        junk = ['ch' num2str(ch_slow(ii))];
    end
    if ~isempty(find(strcmp(var_list,junk)==1))   % Deal with slow channels that are sampled more than once per scan
        eval([junk '= [' junk ' slow_data(:,ii)];'])
    else
        eval([junk ' = slow_data(:,ii);']);
        var_list{n_var} = junk;
        n_var = n_var+1;
    end
end
clear slow_data; 
for ii=1:length(ch_fast)
    jj=find(ch_nums==ch_fast(ii));
    if ~isempty(jj)
        junk = char(ch_names{jj});
    else
        junk = ['ch' num2str(ch_fast(ii))];
    end
    if ~isempty(find(strcmp(var_list,junk)==1))   % Deal with fast channels that are sampled more than once per scan
        eval([junk '= [' junk ' fast_data(:,ii)];']);
    else
        eval([junk '=fast_data(:,ii);']);
        var_list{n_var} = junk;
        n_var = n_var+1;
    end
end
clear fast_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Do some basic conversions: extract the time and date from the %%%%%%
%%%%%% header, and combine split channels                            %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract date/time from the header
time_stamp = datenum(floor(header(:,4:9)));
msec = num2str(header(:,10)/1000,'%0.3f');
msec=msec(:,2:end);
date = datestr(time_stamp,'dd/mm/yyyy');
time = [datestr(time_stamp,'HH:MM:SS') msec];

% Combine split channels (if requested in the setup portion of the routine)
% and slow channels that are sampled more than once per scan
for ii=1:length(var_list)
    if ii>length(var_list),break; end
    varE = var_list{ii};
    if size(eval(varE),2)>1
        junk = eval(varE)';
        junk = junk(:);
        eval([varE '=junk;']);
    end
    if combine_channels && strcmp(varE(end),'E')
        varO = [varE(1:end-1) 'O'];
        if exist(varO,'var')
            Xe = eval(varE); Xo = eval(varO);
            if size(Xo,2)>1
                junk = Xo';
                junk = junk(:);
                Xo=junk;
            end
            jj = find(Xe<0); Xe(jj) = Xe(jj)+B;
            jj = find(Xo<0); Xo(jj) = Xo(jj)+B;
            eval([varE(1:end-1) '= Xo*B + Xe;']);
            var_list{ii} = varE(1:end-1);
            jj = find(strcmp(var_list,varO)==1); 
            if jj==length(var_list), var_list = var_list(1:end-1);
            else var_list = var_list([1:jj-1 jj+1:end]);
            end
            clear(varE,varO,'Xo','Xe');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Save the variables to file. Note that this can be time- %%%%%%
%%%%%% consuming for a large file.                             %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Saving variables to file....');
var_list{length(var_list)+1} = 'header_version';
var_list{length(var_list)+1} = 'header';
var_list{length(var_list)+1} = 'date';
var_list{length(var_list)+1} = 'time';
var_list{length(var_list)+1} = 'fs_slow';
var_list{length(var_list)+1} = 'fs_fast';

if header_version >= 6
    var_list{length(var_list)+1} = 'setupfilestr';
    var_list{length(var_list)+1} = 'sections';
    if isXMP
        var_list{length(var_list)+1} = 'instrument_serial_number_section';
    end
end

eval(['save ' outname '  ' var_list{1} ' -v6 ']); % save the first variable only

for ii=2:length(var_list) % save the rest with append
    eval(['save ' outname '  ' var_list{ii} ' -v6 -append']);
end


if nargout<1, clear var_list; end

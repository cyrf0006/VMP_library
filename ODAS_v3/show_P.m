%% show_P
% Extract and show the average pressure from an ODAS data file.
%%
% <latex>\index{Type B!show\_P}</latex>
%
%%% Syntax
%
%   [P, record] = show_P( fileName )
%
% * [fileName] - String containing the name of the ODAS binary data file.
% * []
% * [P] - Vector of the record-average pressure in physical units.
% * [record] - Vector of the record numbers corresponding to P.
%
%%% Description
%
% Calculate the record-average pressure from a binary ODAS data file. This 
% function can be used to identify which portions of a data file contain data
% that is of interest.  One can then decimate the data file into smaller files 
% that contain those specific portions of data.
%
% For this function to work, the pressure channel address must be present in the
% address matrix.  The configuration file must also contain the appropriate 
% coefficients for converting raw pressure values into physical units.
%
% For legacy ODAS (prior to ODAS v6 data files), the user needs to supply the 
% setup file that was used during data acquisition.  Newer data files can be 
% processed directly if they were obtained using a valid configuration file.
%
% The returned vectors are suitable for plotting the pressure history in a
% data file.
%
%%% Examples
%
%     >> [P records] = show_P( 'my_data_file.p' );
%     >> plot(P); set(gca, 'YDir', 'reverse'); grid on;
%
% Generate a pressure (depth) vs. record (time) line plot.  This plot allows one
% to visually identify the sections of a data file they are interested in.  The
% decimate function can then be used to crop the data file to the portions that
% are of interest.

% *Version History:*
%
% * 2010-01-15 (AWS) support for odas v6 and up
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% * 2012-10-24 (WID) documentation update
% * 2013-02-26 (WID) fix name of pressure channel with setupstr

function [P, record] = show_P(fileName)

% parsing of command line parameters
if nargin < 1; fileName    = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get Input from the User
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(fileName); fileName = getDataFileName; end; % get data file name from the user
if (fileName == 'q'); return;end

% Look for the .p and .mat files, deal with problems
[P,N,E]=fileparts(fileName);
if isempty(E)
   fileName = [fileName '.p'];
elseif (strcmp(E, '.mat'))
    fileName = [N '.p']; % use the binary file
end

% Now find the pressure data (channel 10) in the file
[ch_slow,ch_fast] = query_odas(fileName);

address_matrix = [ch_slow, repmat(ch_fast, size(ch_slow,1), 1)];
address_matrix = address_matrix'; address_matrix = address_matrix(:);
P_index = find(address_matrix == 10);
P_index = P_index(1);
P_step = length(address_matrix);

% We now have all the information needed to extract the pressure data

fid = fopen_odas(fileName,'r'); % safe to assume that file exists.
junk = fread (fid,64,'ushort');
header_size_i    = 18;
record_size_i    = 19;
header_version_i = 11;
setupfile_size_i = 12;
header_size = junk(header_size_i)/2; % in units of 2-byte words
record_size = junk(record_size_i)/2;

data_size = record_size-header_size;

header_version = bitshift(junk(header_version_i), -8) + bitand(junk(header_version_i), 255) /1000;

isXMP = 0;

if (header_version >=6)
    setupfile_size = junk(setupfile_size_i);

    %sanity check to see if there is a valid setup file size in header
    if (setupfile_size <= 0)
        error('incorrect setup file size extracted from first data record');
    end

    setupfilestr = char(fread(fid, setupfile_size))';

    first_record_size = header_size * 2 + setupfile_size;

    instrument_serial_number_section = [];

    cfg = setupstr(setupfilestr);

    if ~isempty(setupstr(cfg, '', 'xmp'))
        isXMP = 1;

        index_to_serial_number = find(address_matrix' == 254)
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

        scan_length = length(address_matrix)
        tmpdbuf = reshape(tmpdbuf,scan_length,data_size/scan_length);

        instrument_serial_number_section = ['xmp_' num2str(tmpdbuf(index_to_serial_number),'%05d')]
    end
else
    first_record_size = record_size*2;
end



fseek(fid, 0, 'eof'); %move to end of file
length_in_bytes = ftell(fid);
total_records = floor((length_in_bytes - first_record_size) / (record_size*2)); % use floor in case their is a rundat the end.


P = zeros(total_records,1); record = P; % pre-assign vectors

fseek(fid, 0, 'bof'); % move to the beginning of the file
fseek(fid, first_record_size, 'bof'); % move to the second record.

for index = (1:total_records)
    fseek(fid, 2*header_size, 'cof'); % skip over the header
    data = fread(fid, data_size, 'short'); 
    P(index) = mean(data(P_index: P_step: end)); % record-average pressure
    record(index) = index;
end

% Now convert into physical units.
if header_version >= 6
    P = convert_odas(P,'P', 'string', setupfilestr, header_version, instrument_serial_number_section);
else
    P = convert_odas(P,'pres','file','setup.txt');
end

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


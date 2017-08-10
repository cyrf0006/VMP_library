%% read_tomi
% Read a range of records from an ODAS data file.
%%
% <latex>\index{Type B!read\_tomi}</latex>
%
%%% Syntax
%
%   [blocks, header, data] = read_tomi( fid, start_block, end_block )
%
% * [fid] file descriptor for an open ODAS data file
% * [start_block] block index of the first segment block
% * [end_block] block index of the last segment block - inclusive
% * []
% * [blocks] vector of the extracted record numbers
% * [header] record header
% * [data] requested record data
%
%%% Description
%
% Read a range of blocks from an ODAS data file. This function is primarily used
% by patch_odas but can also be used directly for the purpose of processing 
% large data files into manageable portions.
% 
%%% Examples
%
%    >> [fid, error] = fopen_odas( 'raw_data_file.p', 'r' );
%    >> if isempty( error ),
%    >>     [blks, header, data] = read_tomi( fid, 1, 1 );
%    >>     fclose( fid );
%    >> end
%
% Read the first block from the data file 'raw_data_file.p'.

% Version History
%
% * 2010-01-15 (AWS) support for odas v6 and up
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-10-23 (WID) documentation added

function [blocks, header, data]=read_tomi(fid,start_block,end_block)

fseek(fid,0,'bof');

junk = fread (fid,64,'ushort');
header_size = junk(18)/2; % in units of 2 byte words
record_size = junk(19)/2;
data_size = record_size - header_size;
header_version_i = 11;

data            = NaN(data_size*(end_block - start_block +1),1);
d_buffer        = ones (data_size,1)*NaN;
header          = ones (header_size*(end_block - start_block +1),1)*NaN;
h_buffer        = ones (header_size,1)*NaN;


header_version = bitshift(junk(header_version_i), -8) + bitand(junk(header_version_i), 255) /1000;

if (header_version >=6)
    setupfile_size = junk(12);
    first_record_size = setupfile_size + header_size * 2;
else
    first_record_size = record_size * 2;
end

fseek(fid, first_record_size + 2*(header_size+data_size)*(start_block-1),'bof');

% there is no data in block zero. 

% Read the file
for index = 1:(end_block - start_block + 1)

% Now read the header
   [h_buffer, count] = fread (fid, header_size,'short');
   if count ~= header_size,
    disp(['Failed to read header at block = ' num2str( start_block + index -1)])
        blocks = index -1;
        header = header(1:(index-1)*header_size);%trim stuff we failed to read
        data   = data  (1:(index-1)*data_size);
        break
   end

% Now read the data
   [d_buffer, count] = fread (fid, data_size,'short');
   if count ~= data_size
    disp(['Failed to read data at block = ' num2str(start_block +index -1)])
        blocks = index -1;
        header = header(1:(index-1)*header_size);%trim stuff we failed to read
        data   = data  (1:(index-1)*data_size);
        break
   end

   data (1+((index-1)*data_size):index*data_size) = d_buffer;      %fill in
   header (1+((index-1)*header_size):index*header_size) = h_buffer;%fill in

blocks = index;
end




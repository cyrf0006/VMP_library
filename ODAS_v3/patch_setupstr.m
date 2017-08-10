%% patch_setupstr
% Patch a configuration string into an existing v6 ODAS data file
%%
% <latex>\index{Type A!patch\_setupstr}</latex>
%
%%% Syntax
%   patch_setupstr( rawDataFile, configFile )
%
% * [rawDataFile] Name of the ODAS data file into which 'configFile' should
%         be embedded.
% * [configFile] Name of the configuration file to embed.
%
%%% Description
%
% Patch an external configuration string into the first record of an ODAS data 
% file.
%
% This function is used to modify the header contained within an ODAS data file.
% One typically performs this task when calibration values contained within a
% data file need to be modified.  It can also be used to convert older data
% files into v6 data files.
%
% When patching a data file, the parameters that directly affect data 
% acquisition *MUST NOT* be changed.  These values include;
%
% # rate,
% # recsize,
% # no-fast,
% # no-slow,
% # matrix.
%
%%% Notes
%
% # The original data file will be backed-up as 'fname_orig.p' only if this file
%   does not already exist.
% # The patched data file will be created as 'fname.p'.
% # A comment will be added to the top of the new setup file string to emphasize
%   that this data file has a patched setup file string. 
% # Certain parameters in the root section must be left unchanged in order to 
%   maintain the original structure of the data. Error messages will be 
%   displayed when attempting to change any of these values.
% # Only the first header will be modified with respect to the new setup file 
%   string size at offset 12.
%
%%% Examples
%
%    >> patch_setupstr( 'data_005.p', 'setup.cfg' );
%
% Embed the configuration string from 'setup.cfg' into the data file 
% 'data_005.p'.  If the backup file 'data_005_orig.p' does not exist, one is 
% created with the contents of 'data_005.p'.
%
%    >> extract_setupstr( 'data_005.p', 'data_005.cfg' );
%    >> edit data_005.cfg
%    >> patch_setupstr( 'data_005.p', 'data_005.cfg' );
%
% Example of how 'patch_setupstr' can be used with 'extract_setupstr' to modify
% the calibration coefficients embedded within a data file.  Hidden in this
% example is how changes are made to the configuration file, 'data_005.cfg', by
% the edit command.
%
%%% Notes
%
% # Older v1 data files can be patched with a configuration file to update them
%   to v6 data files.  This greatly simplifies subsequent data processing and is
%   highly recommended.
% # Patching v1 data files with a v6 header only changes the first header.
%

% *Version History:*
%
% * 2012-03-14 (WID) initial version
% * 2012-03-22 (WID) extensive rewrite, bugs fixed + the following improvements
%                     - error checking when reading, writing, and opening files
% * 2012-03-26 (WID) added support for patching v1 data files with v6 headers
% * 2012-04-11 (WID) replaced inifile_with_instring calls with setupstr
% * 2012-04-25 (WID) updated documentation and added changes for R2007b
% * 2012-09-09 (WID) modified to utilize file_with_ext.  Updated documentation.
% * 2012-11-05 (WID) updated documentation


function patch_setupstr(datafname, setupfname)

[N,M,E,datafname] = file_with_ext( datafname, ...
                       {'' '.p' '.P'}, ...
                       ['Unable to find file: ' datafname] );

[N,M,E,setupfname] = file_with_ext( setupfname, ...
                       {'' '.cfg' '.CFG'}, ...
                       'ODAS v6 configuration file not found.' );

% Attempt to open all required files.  Must find and report all possible 
% errors before proceeding.

% Attempt to open the data file.
[data_fid, msg] = fopen_odas(datafname, 'r');
if (data_fid < 3 || ~isempty(msg))
    error('Unable to open file: %s\n%s', datafname, msg);
end
    
% Attempt to open the setup file.
[setup_fid, msg] = fopen(setupfname, 'r');
if (setup_fid < 3 || ~isempty(msg))
    error('Unable to open file: %s\n%s', datafname, msg);
end

% Read the header, process errors.
[header, size] = fread(data_fid, 64, '*ushort');
if size ~= 64, error(['Unable to read header of file: ' datafname]); end

% set endian format.
if (header(64) == 1), endian = 'l'; else endian = 'b'; end

% set version of the data file.  Should be 0x0001 for legacy data files
if header(11) == 1, version = 1; else version = header(11); end

% Construct new header for v1 data files.
if version == 1, header = convert_header(header); end

% The new generated data file will be created as a temporary file - make 
% sure this temp file doesn't exist - Matlab doesn't guarantee the function call
% 'tempname' returns a unique name.  Loop until one is found.
new_file = tempname;
while exist(new_file, 'file'), new_file = tempname; end
    
[patched_fid, msg] = fopen(new_file, 'w', endian);
if (patched_fid < 3 || ~isempty(msg))
    error('Unable to create file: %s\n%s', new_file, msg);
end
    
% the size in bytes of the setup file string contained in the first record
setupFileSize = header(12);

% read the old and new setup file into single strings
[oldstr, size] = fread(data_fid, double(setupFileSize), '*char*1');
if size ~= setupFileSize, error('Unable to extract setup file.'); end
    
[newstr, size] = fread(setup_fid, inf, '*char*1');
if size == 0, error('Unable to read provided setup file.'); end

oldstr = oldstr';
newstr = newstr';

% skip different safetly checks for v1 data files.  The observed setup file for 
% such files will be invalid - it's the first data block.
if version == 1
    
  % check that the following values are the same:  rate, no-fast, no-slow, 
  % and the matrix
  cfg = setupstr(newstr);
  n_rate = str2double(setupstr(cfg,'root','rate'));
  n_fast = str2double(setupstr(cfg,'root','no-fast'));
  n_slow = str2double(setupstr(cfg,'root','no-slow'));
  n_rows = str2double(setupstr(cfg,'root','num_rows'));
  
  % derive old rate from header...
  freq = header(21) + header(22)/1000;
  o_rate = freq / (header(29) + header(30));
  if abs(n_rate - o_rate) > 0.5
    error('\ndifference in rate detected: old = %d  new = %d', uint16(o_rate), n_rate);
  end
  
  % compare no-fast
  if n_fast ~= header(29)
    error('\ndifference in no-fast detected: old = %d  new = %d', header(29), n_fast);
  end
  
  % compare no-slow
  if n_slow ~= header(30)
    error('\ndifference in no-slow detected: old = %d  new = %d', header(30), n_slow);
  end
  
  % read in the existing matrix.  It is stored as characters in the oldstr 
  % variable.  Remember, endian format matters...
  c = header(29) + header(30);          % c for columns
  r = header(31);                       % r for rows
  oldstr = oldstr(1:c*r*2);             % trim to matrix length (in bytes)
  oldstr = cast(oldstr, 'uint8');       % we want integers, not characters
  oldstr = reshape(oldstr, 2, c*r);     % prepare for conversion into words
  if endian == 'l', msb = 2; lsb = 1; else msb = 1; lsb = 2; end
  matrix = oldstr(msb,:).*256 + oldstr(lsb,:);    % convert bytes to words
  matrix = reshape(matrix, c, r);                 % turn into matrix
  
  % compare the number of rows
  if n_rows ~= r
    error('Invalid number of matrix rows, old: %d new: %d', r, n_rows);
  end
  
  % now read the matrix from the setup file and compare
  rows = setupstr(cfg, 'matrix', 'row[0-9]+');
  
  for r=1:n_rows
    % parse row into integers
    C = textscan(rows{r}, '%d16');

    % compare the number of elements in each row
    if c ~= length(C{1})
      error('Matrix row #%d has wrong number of columns.', r);
    end

    % compare matrix rows.  Note that the matrix is the transpose of what you
    % would expect.
    index = find(C{1} ~= matrix(:,r), 1);
    if ~isempty(index)
      error('Matrix(%d,%d) values differ - old/new = %d/%d', r,index, matrix(index,r), C{1}(index));
    end
  end

else
  % This is a v6 (or larger) data file

  cfg_new = setupstr(newstr);
  cfg_old = setupstr(oldstr);

  n_rate = str2double(setupstr(cfg_new,'root','rate'));
  n_recs = str2double(setupstr(cfg_new,'root','recsize'));
  n_fast = str2double(setupstr(cfg_new,'root','no-fast'));
  n_slow = str2double(setupstr(cfg_new,'root','no-slow'));
  n_rows = str2double(setupstr(cfg_new,'matrix','num_rows'));

  o_rate = str2double(setupstr(cfg_old,'root','rate'));
  o_recs = str2double(setupstr(cfg_old,'root','recsize'));
  o_fast = str2double(setupstr(cfg_old,'root','no-fast'));
  o_slow = str2double(setupstr(cfg_old,'root','no-slow'));
  o_rows = str2double(setupstr(cfg_old,'matrix','num_rows'));

  if n_rate ~= o_rate,
    error('Difference in "rate" detected: old=%d  new=%d\n',o_rate,n_rate);
  end  
  if n_recs ~= o_recs,
    error('Difference in "recsize" detected: old=%d  new=%d\n',o_recs,n_recs);
  end  
  if n_fast ~= o_fast,
    error('Difference in "no-fast" detected: old=%d  new=%d\n',o_fast,n_fast);
  end
  if n_slow ~= o_slow,
    error('Difference in "no-slow" detected: old=%d  new=%d\n',o_slow,n_slow);
  end
  if n_rows ~= o_rows,
    error('Difference in "num_rows" detected: old=%d  new=%d\n',o_rows,n_rows);
  end
    
  % now make sure that the matrix rows are identical 
  n_rows = setupstr(cfg_new,'matrix','row[0-9]+');
  o_rows = setupstr(cfg_old,'matrix','row[0-9]+');

  for i=1:length(n_rows)
    n_row = textscan(n_rows{i}, '%f');
    o_row = textscan(o_rows{i}, '%f');
    if find(n_row{1} ~= o_row{1}),
      error('Matrix differs: #%d, old_row=%s  new_row=%s\n',i,char(o_row),char(n_row));
    end
  end
end     % end safety checks

% show patch date/time at the top of the file to remind user that we are dealing with a patched setup file.
tmpstr = [';' datestr(now,'yyyy-mm-dd HH:MM:SS') ' patched setup file str' char(10) ];
newstr = [tmpstr newstr];
header(12) = length(newstr);

% write the header, return error if not fully written
count = fwrite(patched_fid, header, 'ushort');
if count ~= 64, error('Unable to write header to new data file.'); end
   
% write new setup string, return error if not fully written
count = fwrite(patched_fid, newstr);
if count ~= length(newstr), error('Error writing setup into data file.'); end

ONEMB = 1024*1024;

% It is assumed that the file pointer is positioned at the second record
% before entering this loop.
while ~feof(data_fid)
    [buf, numread] = fread(data_fid,ONEMB);
    fwrite(patched_fid,buf(1:numread));
end

fclose('all');

% Success in creating a patched ODAS data file - time to replace the original
% - rename/move all files to their new locations.
[N,Nfn,M] = fileparts(datafname);
%count = 0;
new_name = [Nfn '_orig.p'];
%while exist(new_name, 'file')
%    count = count + 1;
%    new_name = sprintf('%s%s%d.p', Nfn, '_orig', count);
%end

% Move the old data file.  Don't overwrite the original file if present.
if ~exist(new_name, 'file')
  [success, msg, N] = movefile(datafname, new_name);
  if ~success
      error('Error saving original file to: %s\n%s', new_name, msg);
  end
end

% Move the new data file into the location of the old data file
[success, msg, N] = movefile(new_file, datafname);
if ~success, error('Error moving new file to: %s\n%s', datafname, msg); end

end



function new = convert_header(header)

    % by default, copy all the values.
    new = header;
    
    % data file format to v6.0
    new(11) = bitshift(6, 8);
    
    % setup file size not yet known.  Set to size of data block so we skip over
    % the first empty block when extracting the setup file.
    new(12) = header(19) - header(18);
    
    % product ID is unknown - leave as legacy so we always know this was a 
    % converted file.
    new(13) = 0;
    
    % unknown build number - just leave as 0.
    new(14) = 0;
    
end

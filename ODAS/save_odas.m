%% save_odas
% More efficient way to [re]write a vector in an existing mat file. 
%%
% <latex>\index{Type A!save\_odas}</latex>
%
%%% Syntax
%
%   success_flag = save_odas(file_name, vector_name, vector)
%
% * [file_name] destination .mat file into which the vector should be saved
% * [vector_name] name of the vector being saved
% * [vector] the actual vector that is saved
% * []
% * [success_flag] result of operation, see following table
%
% Possible values for success_flag:
%
% # 0: success, vector found in .mat file is of equal length to the vector 
%      being saved;
% # 1: success, vector not found in .mat file so the vector is saved using the
%      '-append' flag;
% # -1: failure, vector found in .mat file is of different length to the vector
%       being saved.  The .mat file was not modified.
%
%%% Description
%
% This function was written to spare this author recurring frustration
% with the Matlab 'save' command. The creators of Matlab have certainly done
% a marvelous service to the scientific and engineering community. However,
% it appears that they never envisioned that users might want to work with
% large data files and that they might want to update vectors in such a file.
% If you want to append a vector to an existing mat-file, then this
% operation proceeds quite quickly because the vector is simply placed into
% the back of the file. However, if you have loaded a vector, changed its
% value but not its length, and try to update your file by appending this
% modified vector into the mat-file, then Matlab ties itself into a
% multi-dimensional knot. It first copies the entire file into a temporary
% file, empties the original file, copies the vectors and variables
% one-by-one, replaces the target vector when its turn arrives, and then
% continues to copy vectors and variables one-by-one. Microstructure profiles
% tend to have vectors of about 10e6 in size, and re-saving a single vector
% can take up to 2 minutes!. Enter save_odas. This function examines the
% mat-file to see if it already contains a vector with the same name as the one
% that is to be saved. If it does not find a vector of the same name, then it
% simply issues the save  append command.  If it does find a vector with the
% same name, then it checks to see if it has the same length as the one that is
% to be saved. If the lengths are equal, then it uses the binary write command
% to over-write the vector in the mat-file directly. Updating a vector of
% length of about 10e6 takes approx. 1 second. Can you believe it? Can you
% imagine what this has done for my blood pressure?
%
% This function only works for mat-files saved as version 6. You do this by 
% using the v6 flag with the save command. This function, of course does this 
% for you automatically. MathWorks has come out with 2 newer versions of its 
% mat-files, but RSI has no knowledge of the internal structure of these 
% newer-version files.
%
% NOTE: This function works only with real vectors.
%
%%% Examples
%
%    >> vector = 1:1000000;
%    >> save( 'my_big_vector.mat', 'vector' );
%    >> vector(1) = -1;
%    >> save_odas( 'my_big_vector.mat', 'vector', vector );
%
% First creates a .mat file and saves a large vector.  The vector is then 
% modified and the .mat file updated.

% Version History
% * 2000-05-04 (RGL) initial
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-03-21 (RGL) will now save complex vectors by using save directly
%                    with append. Will issue warning rather than error.
% * 2012-10-23 (WID) formatted documentation for publishing

function success_flag = save_odas(file_name,vector_name, vector)

success_flag = -1; % initialize this flag to indicate a failure

if (isreal (vector) ~= 1); 
    warning(['The vector ' vector_name ' must be real']); 
    eval([vector_name ' = vector;']); % copy vector into a variable with the desired name
    eval(['save ' file_name '  ' vector_name ' -append']) % use matlab save command
end % do some input checking
%if (size(vector,1) ~= 1 & size(vector,2) ~= 1); error ('vector must be one-dimensional');end

fid = fopen(file_name,'r+');
if (fid <3); error (['Error opening file = ' file_name]);end
[junk, junk, local_machine_type] = fopen(fid);% Determine native format on this machine

fseek(fid,126,'bof'); % position of the endian flag
endian_flag = fread(fid,1,'int16');
M = 77; % Ascii code for "M"
I = 73; % Ascii code for "I"
if (endian_flag ~= 256*M + I & endian_flag ~= 256*I + M)
   error('Cannot identify machine flag for this M-file')
elseif (endian_flag == 256*M + I); % Local machine type mataches type of file. No further action is needed
elseif (endian_flag == 256*I + M) ; % Local machine type not compatible with file type
   fclose(fid); 
   if (local_machine_type == 'ieee-le'); fid = fopen(file_name,'r+','ieee-be');end % force type match
   if (local_machine_type == 'ieee-be'); fid = fopen(file_name,'r+','ieee-le');end % force type match
end

% find variable name in mat-file. Start by looking for miMATRIX flag and then other parameters
mimatrix = 14; % values assigned my Mathworks
miint8   = 1;
miint32  = 5;
miunit32 = 6;
midouble = 9;
mxdouble_class = 6;

fseek(fid, 0, 'eof'); length_of_file = ftell(fid); % Note where the file ends
fseek(fid,128,'bof'); % position file pointer to the start of the first data element
found_it = -1; end_of_file = -1;

while (found_it == -1 & end_of_file == -1)
   index = ftell(fid); % remember where this data element starts
   type_flag = fread(fid,1,'int32');
   num_of_bytes = fread(fid,1,'int32');
	if (type_flag == mimatrix)% then check for array name and other stuff. Otherwise skip it.
   	if(fread(fid,1,'uint32') == miunit32 & fread(fid,1,'uint32') == 8)
         array_flag = fread(fid, 1, 'uint32'); % read array flag
         if(mod(array_flag, 256) == mxdouble_class); %Check lowest byte only, then go on to read array dimensions
    	     fseek(fid,4,'cof');% skip over undefined region
     	     if (fread(fid,1,'int32') == miint32); % read number of bytes in dimensions array
          	  bytes_in_dimension_array = fread(fid,1,'int32');
               if(bytes_in_dimension_array/4 == 2); % we have an MxN array. Hope that either or both of M and N = 1
                  fseek(fid,bytes_in_dimension_array,'cof');% skip over dimensions
                  test_miint8 = fread(fid,1,'uint32'); %Read array name flag as uncompressed 32-bit word
                  if(test_miint8 > 2^16); % we have data compression
                     fseek(fid, -4, 'cof'); % back up and try again
                     test_miint8 = fread(fid,1,'uint16'); % read as 16-bit word
                     if(test_miint8 == miint8); % we have reached the ascii character string for this data element
                     	length_of_string = fread(fid,1,'uint16'); % the number of characters in the data name string
                        data_name_string = [];
                        data_name_string(1:length_of_string) = fread(fid, length_of_string,'uchar');
                        padding = mod(ftell(fid),8); padding = 8 - padding; padding=mod(padding,8);
                        fseek(fid,padding,'cof'); % skip over padding, if any, for 64-bit boundary
                     end
                  else % we have uncompressed 32-bit data
                     if(test_miint8 == miint8); % we have reached the ascii character string for this data element
                     	length_of_string = fread(fid,1,'uint32'); % the number of characters in the data name string
                        data_name_string = [];
                        data_name_string(1:length_of_string) = fread(fid, length_of_string,'uchar');
                        padding = mod(ftell(fid),8); padding = 8 - padding; padding=mod(padding,8);
                        fseek(fid,padding,'cof'); % skip over padding, if any, for 64-bit boundary
                     end
                  end
                  if(length(vector_name) == length_of_string)
                     if(data_name_string == vector_name) % Eurika, we found it. Pointer is at start of data
                  		found_it = 0;
                        location = ftell(fid);
                     end
                  end
               end
            end
         end
      end
   end
   if (found_it == -1)
      fseek(fid,index,'bof'); % back up to start of this data element
   	fseek(fid,num_of_bytes + 8,'cof'); % skip forward to next data element
      if (length_of_file == ftell(fid)); end_of_file = 0;end % test for end of file
   end
end

if (found_it == 0)
   fseek(fid,location,'bof'); % Move to where variable was found
   test_double = fread(fid,1,'uint32');
   if(test_double == midouble);
      array_length_in_bytes = fread(fid,1,'uint32');
      array_length_in_words = array_length_in_bytes / 8;
      if (length(vector) == array_length_in_words); % we have a size match. Safe to overwrite data
         fseek(fid,0,'cof'); %Dummy move else fwrite does not work
         if(fwrite(fid, vector, 'float64') == length(vector)) % try to insert the data
            success_flag = 0; % indicates successful direct binary write of vector into file
         end
      end
   end
end

if (end_of_file == 0) % then we failed to find the variable
   eval([vector_name ' = vector;']); % copy vector into a variable with the desired name
   eval(['save ' file_name '  ' vector_name ' -append']) % use matlab save command
   success_flag = 1; % indicates successful write using matlab "save" command
end

fclose (fid);

% END of CODE

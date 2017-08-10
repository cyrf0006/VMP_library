%% odas_update
% Perform online update and version check of ODAS Matlab Library
%%
% <latex>\index{Type A!odas\_update}</latex>
%
%%% Syntax
%   [return_str return_num] = odas_update( command )
%
% * [command] Function that should be initiated.  Valid inputs are 'version',
%             'check', and 'update'.
% * []
% * [return_str] Version of the odas_library as a string.  Also returns an empty
%                string with the 'check' command if no update is required.
% * [return_num] Version of the odas_library as a numeric value.
%
%%% Description
%
% This function is used to update the ODAS Matlab Library.  It checks the
% internet to see if the current version of the library is outdated and then
% downloads and installs the newer version if required.
%
% The possible commands are;
%
% # 'version',
% # 'check',
% # 'update';
%
% where 'version' returns a numeric representation of the local version, 
% 'check' goes onto the internet to see if a newer version is available 
% returning the version number if an update is available, or 0 if there is no 
% available update, and 'update' downloads and installs the latest version.
%
%%% Examples
%
%    >> odas_update version
%
% Return the current version of the local ODAS Matlab Library.
%
%    >> if odas_update( 'check' ), 
%    >>     disp('the library should be updated');
%    >> end
%
% Compare the current local version to the online version.  Used to see if an
% update is required.
%
%    >> odas_update update;
%
% Download and install the newest version of the ODAS Matlab Library.  The 
% program will interact with the user to determine the best location for the 
% library.

% Version History:
%
% * 2012-10-26 (WID) initial
% * 2012-12-11 (WID) configured for use on the rocklandscientific.com website.
% * 2012-12-11 (WID) various other fixes, displays options if none given
% * 2013-02-27 (WID) modified for improved version numbers (major,minor,rev).
%                    Broken into two parts - uses the odas_update_remote
%                    function to perform the actuall update.

function [return_str return_num] = odas_update( command )

odas_update_url = 'http://www.rocklandscientific.com/LinkClick.aspx?fileticket=9XvmQL2mJ9A%3d&tabid=55';

version_major =    uint32(3);
version_minor =    uint32(1);
version_revision = uint32(1);

version_odas = uint32( bitor( ...
    bitor( bitshift(version_major, 16), bitshift(version_minor, 8) ), ...
    version_revision ) );

version_odasstr = [num2str(version_major) '.' ...
                   num2str(version_minor) '.' ...
                   num2str(version_revision)];

return_str = '';            % So we don't cause a return error.
return_num = 0;             % So we don't cause a return error.

if nargin == 0,
    disp( ' ' );
    disp( 'Syntax:  odas_update {command} ' );
    disp( ' ' );
    disp( 'Values for {command} include, ' );
    disp( '   version (ver): Display local version number of the ODAS library.' );
    disp( '      check (ck): Check for a newer version online, returns the new version.' );
    disp( '     update (up): Update the ODAS library to the latest version.' );
    disp( ' ' );
    return
end

if strcmpi( command, 'version' ) || strcmpi( command, 'ver' ),
    return_str = version_odasstr;
    return_num = version_odas;
    return
end

if strcmpi( command, 'check' ) || strcmpi( command, 'ck' ),
    
    newstr = check_status;
    
    return_str = '';
    
    if ~isempty(newstr),
        disp('Your version of the ODAS Library is outdated.');
        disp(['Current version = ' version_odasstr]);
        disp(['Online version = ' newstr]);
        return_str = newstr;
    end
    
    return
end

if strcmpi( command, 'update' ) || strcmpi( command, 'up' ),

    status = check_status;
    
    if ~isempty(status),
        remoteFile = acquire_remote_odas_update();
        if isempty(remoteFile), error( 'Download of function failed.' ); end
        [P N E] = fileparts(remoteFile);
        
        % Now add the temp directory to the path
        addpath( P, '-begin' );
        
        savePath = odas_update_remote('update');
     
        rmpath( P );
        rmdir( P, 's' );

        if ~isempty(savePath), savepath; end
    end

    return
end

    function version = check_status()
        % Acquire the remote file, store in a temporary directory
        remoteFile = acquire_remote_odas_update(); 
        if isempty(remoteFile), error( 'Download of function failed.' ); end
        [P N E] = fileparts(remoteFile); 

        % Now add the temp directory to the path
        addpath( P, '-begin' )
        
        [new_str, new_num] = odas_update_remote('version');
        
        rmpath( P );
        rmdir( P, 's' );
        
        version = '';
        
        if new_num > version_odas,
            version = new_str;
        end 
    end


    function download_path = acquire_remote_odas_update()
        % Download the odas_update_remote file and save it in a tmp directory.
        % Return the resulting file - or an empty string if there was an error.
        
        directory = tempname;
        write_status = mkdir( directory );        
        if ~write_status, download_path = ''; return; end
        
        download_path = [directory filesep 'odas_update_remote.m'];
        [download_path,write_status] = urlwrite( odas_update_url, download_path );
        if ~write_status, download_path = ''; end
    end

end
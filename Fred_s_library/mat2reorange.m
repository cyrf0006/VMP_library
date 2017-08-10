
function mat2reorange(profiles, prefix)

% Prepare files for reorange analysis.
% example to run in /home/cyrf0006/WINDEX/data_processing/thorpe-vs-ozmidoz
% > mat2reorange('Tprofiles_riki', 'fine_riki')
%
% in /home/cyrf0006/WINDEX/data_processing/BMix_study:
% > mat2reorange('hit_bottom_tprofiles_slope', 'fine_hitbottom_slope') 

% F. Cyr, Nov. 2011
% --------------------------------------------------------- %

% load profiles*.mat
fid = fopen(profiles);
C = textscan(fid, '%s', 'delimiter', '\n');
list = char(C{1});

no_profiles = size(list, 1);


for i = 1:no_profiles

    fname_in = list(i, :);
    I = find(fname_in==' ');
    fname_in(I)=[];
    fname_out = [prefix '_'];


    if i < 10
        fname_out = [prefix '00' num2str(i) '.dat'];
        disp(fname_out);
    else
        if i < 100
            fname_out = [prefix '0' num2str(i) '.dat'];
                    disp(fname_out);
        else %i>100
            fname_out = [prefix num2str(i) '.dat'];
                    disp(fname_out);
        end
    end
    

    
    load(fname_in)
     
    % -- Computing density and sigma-t for finescale-- %
    DENS = sw_dens(SBS, SBT, P); 
    SIG_T = DENS-1000;
    
    dlmwrite(fname_out, [P, SBT, SBS, SIG_T], 'delimiter',' ','precision',10);
end


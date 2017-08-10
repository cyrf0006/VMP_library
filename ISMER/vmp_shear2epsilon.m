function vmp_shear2epsilon(profiles, zbin)
%usage: vmp_shear2epsilon(no_profiles, 1)
%  or vmp_shear2epsilon('20110922.list', 1)
%
%  Where 'eps_profiles' comes from:
%     "ls -1 profile*.mat | sed 's/\.mat//' > prof_files"
%
%
% modifications:
%
% - F. cyr (2010-11-04)
%    Add a test to ignore treatment when a shear probe is absent
% - F. cyr (2011-01-19)
%    Change "mean" for "nanmean" in the last modif
% - F. Cyr (2011-02-01)
%    Change input method. Input can now be wether the number of
%    profile or a list of profile, same style as with
%    var_profile_cal. Old calling method still work but not
%    supported anymore!
% - F. Cyr (2013-10-22)
%    Added vertical bin size in input (zbin)
% ---------------------------------------------------------- %
tic
% Choose if the function epsilon_vmp.m is explicit or not

answer=0;
while answer==0
    R1 = input('Do you want visual inspection of profiles when suspect? (y/n) ', 's');
    if strcmp('y',R1)==1 
        explicit=1;
        answer=1;
    elseif strcmp('n',R1)==1
        explicit=0;
        answer=1;
    else
        disp('bad answer! please type y or n')
    end
end
        
        
if ischar(profiles)==0
    no_profiles = profiles;
else
    % load profiles*.mat
    fid = fopen(profiles);
    C = textscan(fid, '%s', 'delimiter', '\n');
    pro_files = char(C{1});

    siz = size(pro_files);
    no_profiles = siz(1); %number of profile* files 
end


for count = 1:no_profiles
    
    if ischar(profiles)==0 % the input is a number
        fname_in  = sprintf('profile%3.3i.mat', count);
        fname_out  = sprintf('eps_profile%3.3i', count);
    else % the input is a list.
        fname_in = pro_files(count, :);
        I = find(fname_in==' ');
        fname_in(I)=[];
        fname_out = ['eps_' fname_in];
    end
    
    
    % test to see if a shear probe is absent
    if count==1
        load(fname_in);
        if nanmean(shear1.^2)>0.001 & mean(shear2.^2)>0.001
            disp('both shear probe present')
        elseif nanmean(shear1.^2)>0.001
            shearprobe = 's1';
            disp('S2 absent!?')
        elseif nanmean(shear2.^2)>0.001
            shearprobe = 's2';
            disp('S1 absent!?')
        else
            disp('Both shear probes absent... nothing to do!')
            return
        end
    end
    
    
    d=sprintf('PROCESSING PROFILE %d ...', count);
    disp(d)

    % call shear analysis
    if exist('shearprobe')==0
        vmp_shear_analysis(fname_in,fname_out, zbin, explicit);
    else
        vmp_shear_analysis(fname_in,fname_out, zbin, explicit, shearprobe);
    end
    
end


toc
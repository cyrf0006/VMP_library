function meanTSprofile(ctd_files)

% "ls -1 *.bin > ctd_files"
% 
% For the moment, ctd_files are in 
%   /home/cyrf0006/WINDEX/data_processing
%
% Author: Frederic Cyr, April 2013
% ex:  
    
    
% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

noProfiles = size(ctd,1); 


zbin = 1;
Pbin = [.5:zbin:349.5]';
TMat = nan(length(Pbin), noProfiles);
SMat = nan(length(Pbin), noProfiles);
RMat = nan(length(Pbin), noProfiles);
NMat = nan(length(Pbin), noProfiles);


for i = 1:noProfiles
    fname = ctd(i,:);
    I = find(fname==' ');   
    fname(I) = [];
    disp(fname)
    load(fname)
    
    rho = sw_dens(SBS, SBT, P);
    
    % bin variables
    for ibin = 1:length(Pbin)
        I = find(P >= Pbin(ibin) - zbin/2 & P <= Pbin(ibin) + zbin/2);
        if ~isempty(I) == 1
            TMat(ibin, i) = nanmean(SBT(I));
            SMat(ibin, i) = nanmean(SBS(I));
            RMat(ibin, i) = nanmean(rho(I));
        end
    end
    
    % buoyancy frequency    
 
    N = buoy_freq(rho,P,Pbin,zbin);    
    NMat(:,i) = N;
end

meanT = nanmean(TMat,2);
meanS = nanmean(SMat,2);
meanRho = nanmean(RMat,2);
meanN = nanmean(NMat,2);
outfile = [ctd_files '_out.mat'];
save(outfile, 'Pbin', 'meanT', 'meanS', 'meanRho', 'meanN')
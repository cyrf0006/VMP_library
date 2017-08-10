function boot_NOx_undone(ctd_files, outfile_suffix, zmax, yearmin)

% function boot_NO3_all_stat(ctd_files, outfile_suffix, zmax, yearmin)
%
% THE FUNCTION IS NOT FINISHED SINCE APPARENTLY, TOO FEW BOTTLE ARE
% TAKEN AT DEPTH AND BOOTSTRAP UNDERSAMPLE TOO MUCH THE DATASET. wE
% DECIDED TO PERFOM BIN AVERAGING FIRST.
%
% This function performs bootstrap analysis and monthly average of
% profiles in a list. It also computes the confidence intervals on
% the monthly average. There is 3 main features:
% ex in /home/cyrf0006/PhD/CTD_MPO/Stat16:
%
% >  boot_NOx('dat_profiles','stat16.dat', 500, 1996)
%
% 
%
% "ls -1 *.dat > dat_profiles"
%
% Author: Frederic Cyr, Oct. 2011
% 
% ------------------------------------------------------------- %
global Z S
% -- Some parameters -- %
zbin = 1; % would need to be adjusted if not already binned
zmin=zbin/2;
%zmax=450;
P = [zmin:zbin:zmax]';  
nboot = 500;

% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

n = datenum(ctd(:,1:11)); % in files from ODF2ASCII_bottle.sh


% year restriction
I = find(str2num(datestr(n, 10))<yearmin);
n(I) = [];
bvctd(I,:) = [];

no_files = size(ctd,1); 
NO2_mat = nan(length(P), no_files);
NO3_mat = nan(length(P), no_files);
PO4_mat = nan(length(P), no_files);
zRawVec = [];
NO2RawVec = [];
NO3RawVec = [];
PO4RawVec = [];
for j = 1:length(n)       

    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=5
            disp('Are you sure you have the right file format?')
        end
            
        pres = file(:,1);
        NO2 = file(:,2);
        NO3 = file(:,3);
        PO4 = file(:,4);

        %remove -99
        I = find(pres==-99);
        pres(I) = []; %remove when no pressure
        NO2(I) = [];
        NO3(I) = [];
        PO4(I) = [];
        
        % store raw variables
        zRawVec = [zRawVec; pres];
        NO2RawVec = [NO2RawVec; NO2];
        NO3RawVec = [NO3RawVec; NO3];
        PO4RawVec = [PO4RawVec; PO4];
        
    end
end  % for j


% For raw vectors
I = find(isnan(NO2RawVec)==1);
zNO2 = zRawVec;
zNO2(I) = [];
NO2RawVec(I) = [];
I = find(NO2RawVec==-99);
zNO2(I) = [];
NO2RawVec(I) = [];

I = find(isnan(NO3RawVec)==1);
zNO3 = zRawVec;
zNO3(I) = [];
NO3RawVec(I) = [];
I = find(NO3RawVec==-99);
zNO3(I) = [];
NO3RawVec(I) = [];

I = find(isnan(PO4RawVec)==1);
zPO4 = zRawVec;
zPO4(I) = [];
PO4RawVec(I) = [];
I = find(PO4RawVec==-99);
zPO4(I) = [];
PO4RawVec(I) = [];



% -- bootstrap on tanh fit -- %
disp('bootstrap...')
NO2_boot_b = nan(nboot, 4); % 4 for four coefficients
NO3_boot_b = nan(nboot, 4);
PO4_boot_b = nan(nboot, 4);
NNO2 = length(NO2RawVec); % NO2, NO3, PO4 may not be equal...
NNO3 = length(NO3RawVec);
NPO4 = length(PO4RawVec);
count = 0;
count2 = 0;
for iboot = 1:nboot        
    % NO3 and PO4 only since NO2 has not the same pattern
    if mod(iboot, 100) == 0
        disp(sprintf('iteration no. %d / %d', iboot, nboot))
    end
    r = rand(NNO3,1);
    r = ceil(r*NNO3/1);
    NO3Vec = NO3RawVec(r);       
    [Z, I] = sort(zNO3(r));  
    S = NO3Vec(I);
    a0 = [5 10 .1 100];
    options = optimset('MaxFunEvals', 1500, 'display', 'off');
    [asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
    if exitFlag == 1
        NO3_boot_b(iboot,:) = asol;
    else
        count = count+1;
    end
    
        
    r = rand(NPO4,1);
    r = ceil(r*NPO4/1);
    PO4Vec = PO4RawVec(r);       
    [Z, I] = sort(zPO4(r));  
    S = PO4Vec(I);
    a0 = [1 1 .1 10];
    options = optimset('MaxFunEvals', 1500, 'display', 'off');
    [asol, funVal, exitFlag] = fminsearch ( 'fitFun10', a0, options);
    if exitFlag == 1
        PO4_boot_b(iboot,:) = asol;
    else
        count2 = count2+1;
    end
end

% Result idea
colorscale = 0:.1:1;
figure(1)
clf
count = 1;
plot(NO3RawVec, zNO3, '.k')
set(gca, 'ydir', 'reverse')
hold on
NO3Mat = nan(length(P), nboot);
for i = 1:nboot
    asol = NO3_boot_b(i,:); 
    S3 = asol(1) + asol(2).*tanh((P-asol(3))./asol(4));
    NO3Mat(:,i) = S3';
    
    plot(S3, P, 'color', [1 1 1].*colorscale(count))
    if count == length(colorscale)
        count = 1;
    else
        count = count+1;
    end

end
plot(NO3RawVec, zNO3, '.k')
hold off


disp(['The function is not finished since preliminary results are not satisfying.'])
keyboard



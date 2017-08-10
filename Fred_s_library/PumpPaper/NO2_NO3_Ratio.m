function NO2_NO3_Ratio(ctd_files, outfile_prefix, zbin, zmax, yearmin)
    
% NO2_NO3_Ratio('allStats.list','scatterNO3_all', 10, 500, 2000)
%run in /home/cyrf0006/PhD/Nitrates
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 12 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.02; % top of figure
bots = 0.12; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


% -- Some parameters -- %
%zbin = 10; % would need to be adjusted if not already binned
zmin=zbin/2;
%zmax=450;
P = [zmin:zbin:zmax]';  
nboot = 500;

% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

% for allStats.list
fid = fopen('allStat_date.list');
C = textscan(fid, '%s', 'delimiter', '\n');
date = char(C{1});
n = datenum(date(:,1:11)); % in files from ODF2ASCII_bottle.sh

% for rest
%n = datenum(ctd(:,1:11)); % in files from ODF2ASCII_bottle.sh

% year restriction
I = find(str2num(datestr(n, 10))<yearmin);
n(I) = [];
ctd(I,:) = [];


no_files = size(ctd,1); 
NO2_mat = nan(length(P), no_files);
NO3_mat = nan(length(P), no_files);
PO4_mat = nan(length(P), no_files);
zRawVec = [];
NO2RawVec = [];
NO3RawVec = [];
PO4RawVec = [];
%for j = 1:length(n)       
for j = 1:size(ctd,1)       

    fname = ctd(j, :);
    I = find(fname==' ');   
    fname(I) = [];
    file = load(fname);
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        if size(file, 2)~=5
            disp(['Are you sure you have the right file fortmat? in ' ctd(j,:)])
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
        
        % sort just in case (problem in Stat25, 24-APR-2001_14:44:27.dat)
        [pres, I] = sort(pres);
        NO2 = NO2(I);
        NO3 = NO3(I);
        PO4 = PO4(I);
        
        I = find(diff(pres)==0); %if 2 bottles at same depth
        pres(I) = []; % average them
        NO2(I) = nanmean(NO2([I I+1]));
        NO3(I) = nanmean(NO3([I I+1]));
        PO4(I) = nanmean(PO4([I I+1]));
        NO2(I+1) = []; % remove xtra entry
        NO3(I+1) = [];
        PO4(I+1) = [];

        for i = 1:length(pres)
            [Y, I] = min(abs(pres(i)-P));
            NO2_mat(I,j) = NO2(i);
            NO3_mat(I,j) = NO3(i);
            PO4_mat(I,j) = PO4(i);
        end
    end
end  % for j


% -- Quality control -- %
[X, Z] = meshgrid(1:length(n), P);
zVecNO2 = Z(:);
zVecNO3 = Z(:);
zVecPO4 = Z(:);
NO2Vec = NO2_mat(:);
NO3Vec = NO3_mat(:);
PO4Vec = PO4_mat(:);

I = find(isnan(NO2Vec)==1);
NO2Vec(I) = [];
NO3Vec(I) = [];
zVecNO2(I) = [];
zVecNO3(I) = [];
I = find(NO2Vec==-99);
NO2Vec(I) = [];
NO3Vec(I) = [];
zVecNO2(I) = [];
zVecNO3(I) = [];

I = find(isnan(NO3Vec)==1);
NO3Vec(I) = [];
NO2Vec(I) = [];
zVecNO3(I) = [];
zVecNO2(I) = [];
I = find(NO3Vec==-99);
NO3Vec(I) = [];
NO2Vec(I) = [];
zVecNO3(I) = [];
zVecNO2(I) = [];

keyboard
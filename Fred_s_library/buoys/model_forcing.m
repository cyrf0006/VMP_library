function model_forcing(meteo_files, date_files, t0, tf, varargin);


%function model_forcing(meteo_files, date_files, varargin);
% This function extract and format a climatology to force a
% diffusivity model. The Script is base on the original
% buoy_clim.m, but is now adapted to extract SST and SSS. 
% 
% Function inputs:
%  - meteo_files: a list of meteoXX.dat from 
%     "ls -1 meteo0* | sed 's/\.dat//' > meteo_files" 
%  - date_files: a list of dateXX.dat obtained the same way    
%     "!ls -1 date0* | sed 's/.dat//' > date_files
%  - t0 and tf: initial and final time of the dataset to be
%  created. Must be in mtime format 
%  - varargin: can be whether 'raw', 'hourly', 'daily', 'weekly',
%  'monthly'. The default is 'daily'
%
%  NOTE: meteo*.dat and date*.dat have been extracted from SGDO's
%  raw buoy files from  (ex: MMOB_BOUEE2002_RIMOUSKI_IML4_METOCE.ODF)
%  the method is to 1st use a shell script to modify files from ODF to
%  ASCII: "~/shellscripts/ODF2ASCII.sh  MMOB_BOUEE2008_RIMOUSKI_IML4_METOCE.ODF meteo08 date08"
%
% usage ex: model_forcing('meteo_files', 'date_files', [5 1], [11 30], 'daily');
%    (MAY-1 to NOV-30, daily clim...)
%
% Function outputs: none!
% 
%  The function will create 'daily_forcing.dat' (or raw, weekly,
%  etc.) containing [n   SST   SSS]
%                   ...  ...   ...
%
% author: F. Cyr, feb 2011
%
%
%
% ------------------------------------------------------------------------%

% $$$ % Varargin test
if isempty(varargin)==1
    raw = 0; %save micro = no!
    hourly = 0;
    daily = 1; %default
    weekly = 0;
    monthly = 0;
elseif size(varargin,2)==1
    raw = strcmp(varargin{1}, 'raw');
    hourly = strcmp(varargin{1}, 'hourly');
    daily = strcmp(varargin{1}, 'daily');
    weekly = strcmp(varargin{1}, 'weekly');
    monthly = strcmp(varargin{1}, 'monthly');
else
    disp('Wrong input... try "help model_forcing"')
    return
end


% load meteo_files
fid = fopen(meteo_files);
C = textscan(fid, '%s', 'delimiter', '\n');
meteofiles = char(C{1});

% load date_files
fid = fopen(date_files);
C = textscan(fid, '%s', 'delimiter', '\n');
datefiles = char(C{1});

count = 1;
% store all data in very long vectors (all years in a row)
for i=1:size(meteofiles,1)
    
    system(['rm -rf OUT']);
    system(['sed -s "s/:/ /g" ' datefiles(i,:) ' > OUT']);
    system(['sed -s "s/-/ /g" OUT > /tmp/tmp.txt']);
    system(['mv /tmp/tmp.txt ./OUT']);
    
    meteo = load(meteofiles(i,:));
    datef = load('OUT');
    
    L = length(meteo(:,1)); %number of recorded step
    
    windsp(count:count+L-1) = meteo(:,3);
    windir(count:count+L-1) = meteo(:,7); %degrees
    gustsp(count:count+L-1) = meteo(:,5);
    pres(count:count+L-1) = meteo(:,13); %hPa
    temp(count:count+L-1) = meteo(:,9);
    sst(count:count+L-1) = meteo(:,15); 
    rh(count:count+L-1) = meteo(:,11)./100; % rel hum.
    
    if size(meteo,2)==20 % 2005 and before
        sss(count:count+L-1) = meteo(:,17); 
    elseif size(meteo,2)==22 % after 2005
        sss(count:count+L-1) = meteo(:,19); 
    end
    
    yyyy(count:count+L-1) = datef(:,1);
    mm(count:count+L-1) = datef(:,2);
    dd(count:count+L-1) = datef(:,3);

    clear meteo datef

    count = count+L;
end

%vector with all data time
n = datenum(yyyy, mm, dd);

% Want wind stress?? see buoy_clim.m
% Want heat flux?? see buoy_clim.m

years = unique(yyyy);


%%%%%%%%%%%%%%%%%%
% - daily clim - %
%%%%%%%%%%%%%%%%%%

if daily == 1;
    
    disp('Daily climatology')
    for k = 1:length(years)

        %days in this year
        nn = datenum(years(k), 1, 1):datenum(years(k), 12, 31);
        
        %bissextile (remove feb 29th)
        if length(nn) == 366
            nn(60)=[];
        end

        % restrict to wanted time
        I = find(nn>= datenum(years(k), t0(1), t0(2)) & nn<= datenum(years(k), tf(1), tf(2)));
        nn = nn(I);

        for j = 1:length(nn)
            
            I = find(n == nn(j));

            %daily_wind(j,k) = nanmean(windsp(I)); %windspeed, not tau
            %daily_T(j,k) = nanmean(temp(I)); 
            daily_SST(j,k) = nanmean(sst(I));
            daily_SSS(j,k) = nanmean(sss(I));
            %daily_HR(j,k) = nanmean(rh(I));
            %daily_pres(j,k) = nanmean(pres(I));    
        end
    end

    daily_SST_clim = nanmean(daily_SST, 2);
    daily_SSS_clim = nanmean(daily_SSS, 2); 
    daily_n_clim = [datenum(999, t0(1), t0(2)):datenum(999, tf(1), tf(2))]';

    % save daily clim
    dlmwrite('daily_forcing.dat', [daily_n_clim daily_SST_clim ...
                        daily_SSS_clim],'delimiter',' ','precision',6);
    

end

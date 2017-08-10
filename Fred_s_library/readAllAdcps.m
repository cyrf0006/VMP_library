function readAllAdcps(adcpFiles, gpsFiles, zbin, tbin, zlim,  tlim, outfile)
    
    % readAllAdcps('adcpfiles', 4, 1, [0 100], [datenum(2010, 04, 1) datenum(2012, 10, 31)], 'outifile.mat')
% tbin in minutes!
%
% usage ex: readAllAdcps('adcpfiles_slope', 'gpsfiles_slope', 4, 1, [0 110], [datenum(2010, 04, 1) datenum(2012, 10, 31)], 'ADCPs_1min_4m_slope.mat') 
% readAllAdcps('adcpfiles_slope', 'gpsfiles_slope', 4, 5, [0 110], [datenum(2010, 04, 1) datenum(2012, 10, 31)], 'ADCPs_5min_4m_slope.mat') 
% ** run in /home/cyrf0006/WINDEX/data_processing/workOnAdcps
% Some precautions:
% - in July 2010, we used 2 ADCPs at the same time; thus certainly <------------------ THIS IS VERY IMPORTANT
% have a time vector with repetition if we ask to extract all ADCP
% files... (choose adcFile list consequently)
%
    
% - There are time errors that are corrected in here:
%  * summer 2010: all ADCP files prior to 19 july 2010 (prior 20 July for
%  the one accompaning VMP_G) have 1 hour in advance (i.e. time is
%  UTC+1h). 
%  * 25 october 2012: ADCP is 1min15s late compare to UTC time (t =
%  UTC-1m15s)
%
% modifs: 2013-05-02 - Correct velocities for boat drift. Added
%                      gpsfiles in input. 
%                    - same magnetic variation for all files


mag_var = -18.5;
    
zVec = zlim(1)+zbin/2:zbin:zlim(2)-zbin/2;
tbin = tbin/(24*60);
timeVec = [];
adcpDepth = 1; % We assume depth = 1m
%timeVec = tlim(1):tbin:tlim(2);

Ugrid = [];%sparse(length(zVec), length(timeVec));
Vgrid = [];%sparse(length(zVec), length(timeVec));


fid = fopen(adcpFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
adcpFiles = char(C{1});
noFiles = size(adcpFiles, 1);

fid = fopen(gpsFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
gpsFiles = char(C{1});

for i = 1:noFiles
    
    fname = adcpFiles(i, :);
    I = find(fname==' ');   
    fname(I) = [];    
    
    %disp([fname])
    
    ADCP = rdradcp(fname, 1);
    %load(fname)
    

    U = ADCP.east_vel;
    V = ADCP.north_vel;
    mtime = ADCP.mtime;
    z = ADCP.config.ranges+adcpDepth;
         
    % remove last timestep (often reading errors)
    U(:,end) = [];
    V(:,end) = [];
    mtime(end) = [];
    
    % Time correction (see info in help menu)
    if mtime(1) > datenum(2010, 01, 01) & mtime(1) < datenum(2010, 07, 20);
        disp('time correction t = t-1h')
        mtime = mtime - 1/24;
    end
    
    if strcmp(fname, ['../mission_juil2010/VMP_G/2010_07_19/ADCP/' ...
                      'D300_000.000']) == 1
        disp('time correction t = t-1h')
        mtime = mtime - 1/24;   
    end
        
    if mtime(1) > datenum(2012, 10, 24) & mtime(1) < datenum(2012, 10, 26);
        disp('time correction t = t+1m15s')
        mtime = mtime + 1.25/24/60;
    end   
    
    
    ti = round(datenum(mtime(1)*24*60))/24/60;
    tf = round(datenum(mtime(end)*24*60))/24/60;     
    timeV = ti:tbin:tf;
    Utmp = nan(size(U,1), length(timeV));
    Vtmp = nan(size(U,1), length(timeV));
    Utmp2 = nan(length(zVec), length(timeV));
    Vtmp2 = nan(length(zVec), length(timeV));    
    
    
    for j = 1:length(timeV)
        I = find(mtime>=timeV(j)-tbin/2 & mtime<=timeV(j)+tbin/2);
        Utmp(:,j) = nanmean(U(:,I), 2);
        Vtmp(:,j) = nanmean(V(:,I), 2);
    end
      
    
    for j = 1:length(zVec)
        I = find(z>=zVec(j)-zbin/2 & z<=zVec(j)+zbin/2);
        Utmp2(j, :) = nanmean(Utmp(I, :), 1);
        Vtmp2(j, :) = nanmean(Vtmp(I, :), 1);
    end
    
    % ---- last corrections ---- %
    % Correct for magvar
    if ADCP.config.magnetic_var < -20 | ADCP.config.magnetic_var > -18
        [Utmp2, Vtmp2] = rotate_vecd(Utmp2, Vtmp2, -mag_var);
    end
    
    % correct for boat drift
    gpsFname = gpsFiles(i, :);
    I = find(gpsFname==' ');   
    gpsFname(I) = [];
    
    Utmp5 = Utmp2;   
    [Utmp2, Vtmp2] = adcp_velcor(Utmp2, Vtmp2, timeV, gpsFname);
% $$$     figure(1)
% $$$     clf
% $$$     imagesc(Utmp2-Utmp5);
% $$$     caxis([-1 1])
% $$$     colorbar
% $$$     title(datestr(timeV(1)))
% $$$     keyboard
    % ------------------------------- %
    

     % store data
    Ugrid = [Ugrid Utmp2];
    Vgrid = [Vgrid Vtmp2];
    timeVec = [timeVec timeV];
end

% Remove NaNs-columns
I = nansum(Ugrid, 1);
J = find(~isnan(I)==1);
Ugrid = Ugrid(:,J);
Vgrid = Vgrid(:,J);
timeVec = timeVec(J);


save(outfile,'Ugrid', 'Vgrid', 'timeVec', 'zVec');    
    
    
    
    
    
function TSG2map(TSG_FILES, latLims, lonLims, months)
    
% usage ex: TSG2map('tsg.list', [45 52], [-72 -55], [1:12])



% Preambule:  Define regular vectors/grid
dlat = .01; % ~<1km
dlon = .02; % ~<1km
latVec = latLims(1):dlat:latLims(end);
lonVec = lonLims(1):dlon:lonLims(end);
    
[lonGrid latGrid] = meshgrid(lonVec, latVec);

    
fid = fopen(TSG_FILES);
C = textscan(fid, '%s', 'delimiter', '\n');
tsgFiles = char(C{1});

noFiles = size(tsgFiles,1);


timeStored = [];
latStored = [];
lonStored = [];
tempStored = [];
salStored = [];


% Open and store all files
for i = 1:noFiles

    fname = tsgFiles(i, :);
    disp(fname)
    I = find(fname==' ');   
    fname(I) = [];
       
    fid = fopen(fname);
    C = textscan(fid, '%s %f %f %f %f %f');
    fclose(fid);

    n = datenum(char(C{1}), 'dd-mm-yyyy_HH:MM:SS.00');
    %n = datenum(char(C{1}), 'dd-mmm-yyyy_HH:MM:SS.00');


% $$$     % Bin time in 1-minute interval
% $$$     %timeVec = min(n):60*5/86400:max(n);
% $$$     if n(2)-n(1)<59/86400
% $$$         disp('short interval!!')
% $$$         keyboard
% $$$     end
    
    timeStored = [timeStored; n];
    latStored = [latStored; C{2}];
    lonStored = [lonStored; C{3}];
    tempStored = [tempStored; C{4}];
    salStored = [salStored; C{5}];
end
disp(sprintf(' -> %d files read', noFiles))


disp('Quality tests')
% Quality test (bad data - 17-Nov-1858)
I = find(str2num(datestr(timeStored, 10))<2000);
timeStored(I) = [];
latStored(I) = [];
lonStored(I) = [];
tempStored(I) = [];
salStored(I) = [];

I = find(tempStored<-2);
timeStored(I) = [];
latStored(I) = [];
lonStored(I) = [];
tempStored(I) = [];
salStored(I) = [];

I = find(salStored>40);
timeStored(I) = [];
latStored(I) = [];
lonStored(I) = [];
tempStored(I) = [];
salStored(I) = [];
disp(' -> passed!')


% monthly means
disp('Monthly climatologies')
for i = 1:length(months) 
    outfile = ['TSGclimato_' datestr(datenum(1,months(i),1),5)];
    disp([' -> proceeding to ' outfile '.mat']);
    
    I = find(str2num(datestr(timeStored, 5))==months(i));
    if isempty(I) == 1
        disp(['No data to process in ' datestr(datenum(1,months(i),1),3)])
    else
        lat = latStored(I);
        lon = lonStored(I);
        T = tempStored(I);
        S = salStored(I);
        
        TMat = nan(size(latGrid));
        SMat = nan(size(latGrid));
        for j = 1:length(latVec)
            for k = 1:length(lonVec)
                J = find(lat>=latVec(j)-dlat/2 & lat<latVec(j)+dlon/2 ...
                         & lon>=lonVec(k)-dlon/2 & lon<lonVec(k)+dlat/2); 
                if isempty(J)~=1
                    TMat(j,k) = nanmean(T(J));
                    SMat(j,k) = nanmean(S(J));
                end
            end
        end

        save(outfile, 'TMat', 'SMat', 'latVec', 'lonVec')      
        
    end
end



% $$$ keyboard
% $$$ 
% $$$ 
% $$$ figure(2)
% $$$ clf
% $$$ m_proj('mercator', 'long',[min(lonVec) max(lonVec)],'lat',[min(latVec) max(latVec)]);
% $$$ m_grid('box','fancy')
% $$$ %m_contourf(lonVec, latVec,TMat, [-2:.1:20], 'linestyle', 'none')     
% $$$ m_pcolor(lonVec, latVec,TMat)     
% $$$ shading interp
% $$$ 
% $$$ hold on
% $$$ m_gshhs_h('patch',[.5 .5 .5]); %coastlines
% $$$ hold off

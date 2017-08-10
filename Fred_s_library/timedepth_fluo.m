function timedepth_fluo(profileList, varargin)

% function timedepth_fluo(profileList, varargin)
%
% varargin is respectively zmax and zbin
%
% ex: timedepth_eps('20110714.list', 100, .2)
% where:
% "ls profile_2011-07-14_0* > 20110714.list"
%
% 
%

% Varargin test
if isempty(varargin)==1
    zmax = 324.5;
    zbin = 1;
    zmin = 0.5;
elseif size(varargin,2)==1
    zmax = varargin{1};
    zbin = 1;
    zmin = zbin./2;
elseif size(varargin,2)==2
    zmax = varargin{1};
    zbin = varargin{2};
    zmin = zbin./2;
else
    disp('Wrong input... try "help var_profile_cal"')
    return
end


% ------- for water column stability -------- %
ddz = .25;
P_N2 = ddz/2:ddz:max(zmax);
g = 9.81;
% ------------------------------------------- %


% load eps_files file
fid = fopen(profileList);
C = textscan(fid, '%s', 'delimiter', '\n');
prof_files = char(C{1});

no_profile = size(prof_files, 1); %number of *.P files 


%zbin=1;
zvec = [zmin:zbin:zmax]';
rhoMat = nan(size(zvec,1), no_profile);
TMat = nan(size(zvec,1), no_profile);
fluoMat = nan(size(zvec,1), no_profile);


%% loop on profiles %%
for iprofile = 1:no_profile

    fname = prof_files(iprofile, :);
    I = find(fname==' '); %remove white space  
    fname(I) = [];
    disp(fname);
    
    load(fname);
    
    % time of survey
    if iprofile==1
        t0 = mtime(1);
    elseif iprofile==no_profile
        tf = mtime(end);
    end
    

    % bin and store slow data
    rho = sw_dens(SBS, SBT, P);
    for i = 1:length(zvec)
        I = find(P >= zvec(i)-zbin/2 & P < zvec(i)+zbin/2);
        TMat(i,iprofile) = nanmean(SBT(I));
        rhoMat(i,iprofile) = nanmean(rho(I));
    end

    % bin and store fast data
    for i = 1:length(zvec)
        I = find(p >= zvec(i)-zbin/2 & p < zvec(i)+zbin/2);
        fluoMat(i,iprofile) = nanmean(fluoro(I));
    end    
        
    timeVec(iprofile) = mtime(1);
    
end 



% Fluo only
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
contourf(timeVec, zvec, fluoMat, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
cb = colorbar;
ylabel(cb,'fluorescence (ppb', 'FontSize', 10)
datetick('x', 15)
title(datestr(timeVec(1),1))
ylabel('depth (m)')
cax = get(gca, 'clim');

% Fluo + isopycnals
figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 35 15])
contourf(timeVec, zvec, fluoMat, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
cb = colorbar;
ylabel(cb,'fluorescence (ppb', 'FontSize', 10)
datetick('x', 15)
title(datestr(timeVec(1),1))
ylabel('depth (m)')
hold on
V = 1018:0.2:1028;
contour(timeVec, zvec, rhoMat, V,'k')
caxis(cax)

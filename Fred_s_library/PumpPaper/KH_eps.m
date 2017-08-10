clear

load ~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/ADCP/adcp_2009-10-01.mat;

load ~/WINDEX/data_processing/Mission_tadoussac_2009/eps_plot_info.mat;


t0 = datenum(2009, 10, 01, 17, 30, 0);
tf = datenum(2009, 10, 01, 19, 15, 0);
%t0 = datenum(2009, 10, 01, 17, 47, 0);
%tf = datenum(2009, 10, 01, 18, 40, 0);
xticks = t0:15/60/24:tf;



intens = squeeze(ADCP.intens(:,4,:));
zVec = ADCP.config.ranges;
timeVec = ADCP.mtime;
zmax = max(zVec);

I = find(timeVec>=t0 & timeVec <=tf);
timeVec = timeVec(I);
intens = intens(:,I);

% find bottom
I = find(zVec>15);
[Y, J] = max(intens(I,:));
J = J+(size(intens,1)-length(I));
I = find(Y<120);
J(I) = length(zVec);
bot = nan(1, length(J));
for i = 1:length(J)
    bot(i) = zVec(J(i));
end
bot = runmean(bot, 3);

% echo attenuation correction
for i = 1:length(timeVec)
    intens(:,i) = intens(:,i).*20.*log10(zVec);
end


I = find(timevec>=t0 & timevec <=tf);
eps_mat = eps_mat(:,I);
K_mat = K_mat(:,I);
rho_mat = rho_mat(:,I);
timevec = timevec(I);


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[10 10 25 5])

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.05 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


imagesc(timeVec, zVec, intens)
colormap(flipud(bone))

%caxis([50 150])
caxis([1000 3000])
freezeColors
drawnow
hold on


% dissipation
field=K_mat;
x=timevec;
y=zvec;
epsMin = 1e-7; % for eps.
epsMax = 1e-4;
epsMin = 1e-5; % for K
epsMax = 1e-1;
rectbar = [0.91 0.15 0.01 0.75]; 
for j = 1:5:size(field,2)     
    data = field(:,j);
    x(1:length(data)) = timevec(j);
    y = zvec;      
    if j == 1 % with cbar
        cb = colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),1,rectbar, x'.*0);
    else % nocbar
        colour_profile(x',y',log10(data),log10(epsMin),log10(epsMax),0,rectbar, x'.*0);
    end 
end


% isopycnals
%contour(timevec, zvec, rho_mat, 15, 'color', [1 1 1], 'linewidth', 1)
set(gca, 'xtick', xticks)
xlim([t0 tf])
%datetick('x', 15, 'keeplimits')
set(gca, 'xticklabel', datestr(xticks, 15))
ylabel(cb, 'K (m^2 s^{-1})')

% bottom
x = [timeVec(1) timeVec  timeVec(end) timeVec(1)];
y = [zmax bot zmax zmax];
patch(x, y, [1 1 1]*.6)
hold off

ylabel('Depth(m)')
adjust_space

set(gcf, 'renderer', 'painters')
print('-depsc2', 'KH_eps.eps')

clear

load ~/WINDEX/data_processing/Mission_tadoussac_2009/2009_10_01/ADCP/adcp_2009-10-01.mat;

% $$$ t0 = datenum(2009, 10, 01, 17, 40, 0);
% $$$ tf = datenum(2009, 10, 01, 18, 50, 0);
% $$$ t0 = datenum(2009, 10, 01, 17, 47, 0);
% $$$ tf = datenum(2009, 10, 01, 18, 40, 0);




intens = squeeze(ADCP.intens(:,4,:));
zVec = ADCP.config.ranges;
timeVec = ADCP.mtime;

% $$$ I = find(timeVec>=t0 & timeVec <=tf);
% $$$ timeVec = timeVec(I);
% $$$ intens = intens(:,I);

% find bottom
I = find(zVec>15);
[Y, J] = max(intens(I,:));
J = J+(size(intens,1)-length(I));
I = find(Y<130);
J(I) = length(zVec);
bot = nan(1, length(J));
for i = 1:length(J)
    bot(i) = zVec(J(i));
end

bot = runmean(bot, 5);
J = find(abs(bot-max(zVec))<.01);
bot(J) = []; % no topography means cannot see the bottom
timeVec(J) = []; % no topography means cannot see the bottom

save HLC_topo_2009-10-01.mat timeVec bot
function boot_o2_all(ctd_files)

% function boot_ctd(ctd_files)
%
% This function performs bootstrap analysis and monthly average of
% profiles in a list. It also computes the confidence intervals on
% the monthly average. There is 3 main features:
% 

% "ls -1 *.bin > ctd_files"
%
% Author: Frederic Cyr, Oct. 2011
% 
% ------------------------------------------------------------- %

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 15 15])

% -- Some parameters -- %
zbin = 1; % would need to be adjusted if not already binned
zmin=1;
zmax=330;
P = [zmin:zbin:zmax]';   
nboot = 500;

% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

no_files = size(ctd,1); 

% compute time vector for considered profiles
% Only work for Peter's Perl renaming script
n = datenum(str2num(ctd(:,1:4)), str2num(ctd(:,6:7)), str2num(ctd(:,9:10)),str2num(ctd(:,12:13)), str2num(ctd(:,15:16)), 0 );

for j = 1:length(n)       
    
    file = load(ctd(j,:));
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        if size(file, 2)~=5
            continue
        end
        
        pres = file(:,1);
        oxyg = file(:,5);
        
        %remove -99
        I = find(oxyg==-99);
        oxyg(I) = NaN;

        % bin profiles
        for k = 1:length(P)
            I = find(pres>(zmin+(k-1)*zbin) & pres<(zmin+(k)*zbin));
            oxyg_bin(k, :) = nanmean(oxyg(I));
        end
        oxyg = oxyg_bin;
        pres = [round(min(pres)):zbin:round(max(pres))]';
        
        if max(pres) >= max(P) & min(pres)<=min(P)
            I = find(pres>=min(P) & pres<=max(P));
            
            % fill matrix with empty value
            O2(P,j)=NaN;
            
            O2(pres(I),j) = oxyg(pres(I));

        else % CTD profile doesnt start at p=1m or pmax smaller than 300m
            
            %I = find(file(:,1)>=max(min(pres), min(P)) & file(:,1)<=min(max(pres), max(P)));
            I = find(pres>=max(min(pres), min(P)) & pres<=min(max(pres), max(P)));
            % fill matrix with empty value
            O2(P,j)=NaN;
            
            % substituate empty value where it is possible
            O2(pres(I),j) = oxyg(pres(I));
        end
        
    else %file empty
        O2(P,j)=NaN;
    end
    
end  % for j

% -- Quality control -- %
% Remove empty columns
I = find(sum(~isnan(O2),1)>0);
O2 = O2(:,I);

% ml/l to mg/l
O2 = O2.*1.42903;
% mg/l t umol/l
O2 = O2.*1000/32;

% Eliminate P > 318m
I = find(P>318);
O2(I,:) = nan(zmax-318, size(O2,2));

% vertical interpolation to remove NaNs
for j = 1:size(O2,2)
    I = find(isnan(O2(:,j))==1);
    if ~isempty(I)==1 
        p = P;
        o = O2(:,j);
        o(I) = [];
        p(I) = [];
        O2(:,j) = interp1(p, o, P);
    end
end

% -- bootstrap -- %
disp('bootstrap...')
NO = size(O2,2); % NT and NS may not be equal...
for b = 1:nboot        
    r = rand(NO,1);
    r = ceil(r*NO/1);
    O2_boot_b(:,b) = nanmean(O2(:,r),2);
end

% --- mean values --- %
O2_boot_dot = nanmean(O2_boot_b,2);


% --- Compute 95% confidence interval --- %
O2_sort = sort(O2_boot_b, 2);

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

O2_2p5 = O2_sort(:,CI_2p5);
O2_97p5 = O2_sort(:,CI_97p5);

O2_ave = O2_boot_dot;




% Oxygen plot
figure(1)
plot(O2_ave, P, 'k', 'linewidth', 1)
hold on

plot([62.5 62.5], [150 350], '--k')

x1 = O2_97p5;
x2 = O2_2p5;

nann = find(isnan(x1)==1);
if ~isempty(nann)==1 % remove NaNs
    
    x1(nann)=[];
    x2(nann)=[];
    PP = P;
    PP(nann)=[];
    patch([x1; flipud(x2); x1(1)], [PP; flipud(PP); PP(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(O2_ave, P, 'k', 'linewidth', 1)
    hold off
else
    patch([x1; flipud(x2); x1(1)], [P; flipud(P); P(1)], [.8 .8 .8], 'edgecolor', 'none');
    plot(O2_ave, P, 'k', 'linewidth', 1)
    hold off
end


set(gca, 'ydir', 'reverse')
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
set(gca, 'yminortick', 'on')
set(gca, 'xminortick', 'on')
set(gca, 'box', 'on')
axis([50 350 0 320])
ylabel('P (dbar)')    
xlabel('diss. O_2(\mu mol/l)')

dlmwrite('boot_O2_all.dat', [P O2_ave O2_2p5 O2_97p5],'delimiter',' ','precision',6)



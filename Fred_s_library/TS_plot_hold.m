function TS_plot_hold(fnames)

% TS_plot_hold(fnames)
% usage ex: TS_plot_hold(['profile035';'profile040'; 'profile060'; 'profile065'])
colorlegend = ['b'; 'r'; 'k'; 'y'; 'm'; 'g'];
%colorlegend = ['b'; 'r'; 'm'; 'g'];
file =1;
%keyboard

S = size(fnames);
sub=1;

%%%%%%%%%%%%%%%
%   SUBPLOT   %
%%%%%%%%%%%%%%%

% while file<=S(1)
%     
% load(fnames(file, :))
% 
% dat = date(1,:);
% h = str2num(time(1, 1:2));
% m = str2num(time(1, 4:5)); 
% s = str2num(time(1, 7:8)); 
% dd = str2num(dat(1:2));
% mm = str2num(dat(4:5));
% yyyy = str2num(dat(7:10));
% 
% n = datenum(yyyy,mm,dd,h,m,s);
% 
% L(file, :) = datestr(n,15);
% 
% subplot(1,S(1),sub)
% plot(SBS,SBT, colorlegend(file))
% title({'Diagramme T-S'; datestr(n,0)}, 'FontSize', 15)
% xlabel('Salinite (psu)', 'FontSize', 15)
% ylabel('T (degC)', 'FontSize', 15)
% axis([26 35 0 6])
% set(gca, 'FontSize', 15)
% 
% sub=sub+1;
% file = file+1;
% 
% end


%%%%%%%%%%%%%%%
%  HOLD PLOT  %
%%%%%%%%%%%%%%%

while file<=S(1)
    
load(fnames(file, :))

dat = date(1,:);
h = str2num(time(1, 1:2));
m = str2num(time(1, 4:5)); 
s = str2num(time(1, 7:8)); 
dd = str2num(dat(1:2));
mm = str2num(dat(4:5));
yyyy = str2num(dat(7:10));

n = datenum(yyyy,mm,dd,h,m,s);

%L(file, :) = datestr(n,15); %1-day plot
L(file, :) = datestr(n,1); %multiple-day plot

plot(SBS,SBT, colorlegend(file))
%title({'Diagramme T-S'; datestr(n,1)}, 'FontSize', 15) %1-day plot
%title('Diagramme T-S', 'FontSize', 14) %multiple-day plot
xlabel('Salinite (psu)', 'FontSize', 14)
ylabel('T (^{\circ}C)', 'FontSize', 14)
axis([26 34 -1 9])
set(gca, 'FontSize', 14)

sub=sub+1;
file = file+1;
hold on
end


legend(L, 'Location','SouthWest')
hold off


function odf2mat_buoy(infile, outfile)

 %!rm -f databuoy datebuoy
% spaces are important!    
system(['sh ~/shellscripts/velbuoy_ODF2ASCII.sh ' infile ' databuoy datebuoy']);

adcp_raw = load('databuoy'); % z, u, u_f, v, v_f, w, w_f, error


% this part only work if CTD  filename convention is respected
fid = fopen('datebuoy');
C = textscan(fid, '%s', 'delimiter', '\n');
date_str = char(C{1}); 

!rm databuoy datebuoy

% isolate depth vector
count = 1;
while adcp_raw(count+1, 1)-adcp_raw(count, 1)>0
    z(count) = adcp_raw(count, 1);
    count = count+1;
end
z(count) = adcp_raw(count, 1);
z = round(z)';

date_str = date_str(1:length(z):length(date_str), :);
mtime = datenum(date_str);


east_vel = adcp_raw(:,2);
east_vel_error = adcp_raw(:,3);
I = find(east_vel_error==3);
east_vel(I)=NaN;

north_vel = adcp_raw(:,4);
north_vel_error = adcp_raw(:,5);
I = find(north_vel_error==3);
north_vel(I)=NaN;

vert_vel = adcp_raw(:,6);
vert_vel_error = adcp_raw(:,7);
I = find(vert_vel_error==3);
vert_vel(I)=NaN;

east_vel = reshape(east_vel, length(z), length(east_vel)/length(z));
north_vel = reshape(north_vel, length(z), length(north_vel)/length(z));
vert_vel = reshape(vert_vel, length(z), length(vert_vel)/length(z));


% here we should make a struct... 
% But for now it is fine...


%dlmwrite(outfile, vert_vel, 'delimiter',' ','precision',6);


save(outfile, 'mtime', 'z', 'east_vel', 'north_vel', 'vert_vel')





function [T, S, n] = ctd_matrix(ctd_files)
% build a matrix with all profiles in a list

disp('Extracting all CTD profiles')
    
% -- Some parameters -- %
dz = 1; % would need to be adjusted if not already binned
P = [1:dz:300]';   


% -- Getting infos on profiles -- %
% load file names
fid = fopen(ctd_files);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});

no_files = size(ctd,1); 

% compute time vector for considered profiles
% Only work for Peter's Perl renaming script
n = datenum(str2num(ctd(:,1:4)), str2num(ctd(:,6:7)), str2num(ctd(:,9:10)),str2num(ctd(:,12:13)), str2num(ctd(:,15:16)), 0 );


T = nan(length(P), no_files);
S = nan(length(P), no_files);

for i = 1:no_files
        
    fname=ctd(i,:);
    I = find(fname==' ');   
    fname(I) = [];  
    file = load(fname);
    
    % Quality control + save into matrix T and S
    if ~isempty(file)==1 % check if file empty
        
        pres = file(:,1);
        temp = file(:,2);
        sali = file(:,3);
        
        %remove -99
        I = find(temp==-99);
        temp(I) = NaN;
        I = find(sali==-99);
        sali(I) = NaN;
        
        if max(pres) >= max(P) & min(pres)==min(P)
            I = find(file(:,1)>=min(P) & file(:,1)<=max(P));
            
            % fill matrix with empty value
            T(P,i)=NaN;
            S(P,i)=NaN;
            
            T(pres(I),i) = temp(I);
            S(pres(I),i) = sali(I);               

        else % CTD profile doesnt start at p=1m or pmax smaller than 300m
            
            I = find(file(:,1)>=max(min(pres), min(P)) & file(:,1)<=min(max(pres), max(P)));
            % fill matrix with empty value
            T(P,i)=NaN;
            S(P,i)=NaN;
            
            % substituate empty value where it is possible
            T(pres(I),i) = temp(I);
            S(pres(I),i) = sali(I);
            
        end
        
    else %file empty
        T(P,i)=NaN;
        S(P,i)=NaN;
    end
    
end  % for i


% -- Quality control -- %
% Remove empty columns 
% (remove both T and S if one column is empty)
I = find(sum(~isnan(T),1)>0);
T = T(:,I);
S = S(:,I);
n = n(I);
I = find(sum(~isnan(S),1)>0);
S = S(:,I);
T = T(:,I);
n = n(I);


% vertical interpolation to remove NaNs
for j = 1:size(T,2)
    I = find(isnan(T(:,j))==1);
    if ~isempty(I)==1 
        p = P;
        t = T(:,j);
        t(I) = [];
        p(I) = [];
        T(:,j) = interp1(p, t, P);
    end
end
for j = 1:size(S,2)
    I = find(isnan(S(:,j))==1);
    if ~isempty(I)==1
        p = P;
        s = S(:,j);
        s(I) = [];
        p(I) = [];
        S(:,j) = interp1(p, s, P);
    end
end


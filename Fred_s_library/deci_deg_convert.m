function New_coord = deci_deg_convert(Coordfile)

% function New_coord = deci_deg_convert(Coord)
disp('signs not working...')
Coord=load(Coordfile);
    
rest = abs(Coord)-abs(floor(Coord));

for i = 1:size(rest,2)
    I = find(rest(:,i)<0);
    if ~isempty(I)==1
        rest(I,i)=1+rest(I,i);
    end
end


min = rest.*60;




New_coord = [floor(Coord(:,1)) min(:,1) floor(Coord(:,2)) min(:,2)];


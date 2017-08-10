function is = isthere(lat1, lon1, lat2, lon2, rad)

%  is = isthere(lat1, lon1, lat2, lon2, rad)
%
% Tells if the distance between 2 pts is within a certain radius
% "rad". Thus if the distance between [lat1, lon1] and [lat2, lon2]
% is larger than rad
%
% Returns is = 1 if the point is considered there, is = 0
% otherwise!
%
% Frederic Cyr (2011-02-07)
% -------------------------------------------------------- %

is = 0;
range = m_lldist([lon1, lon2], [lat1, lat2]);

if range <= rad
    is=1;
end







    

  


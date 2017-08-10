function LineHndls = lines(varargin)
%
% lines.m--Like line.m, but line is coloured according to value of c...
% c should be scaled to some colourmap (i.e., red to blue? jet?)
% Tricky part is getting colourmap not to get overloaded...use 
% get(0,factoryFigureDithermap)?
%
% Syntax: LineHndls = lines(x,y,<z>)
%
% e.g.,   t = -pi:pi/500:pi; x=sin(5*t); y=cos(3*t); z=t; 
%         %c=log([1:length(t)]);
%         c=(t+pi).^0.5;
%         LineHndls = lines(x,y,z,c,spring(64));

% Developed in Matlab 6.1.0.450 (R12.1) on SUN OS 5.8.
% Kevin Bartlett(bartlett@soest.hawaii.edu), 2002/02/22, 13:38
%------------------------------------------------------------------------------

% Parse the input arguments.
if nargin>5,
   error([mfilename '.m--Too many input arguments.']);
end % if 

x = varargin{1};
y = varargin{2};

for ArgCount = 3:nargin,

   if all(size(varargin{ArgCount}(:))==size(x(:))),

      if exist('z')~=1,
         z=varargin{ArgCount};
      elseif exist('c')~=1,
         c=varargin{ArgCount};
      else
         error([mfilename '.m--Line colour specification must be nx3']);
      end % if

   else

      if exist('CMap')==0,
         CMap = varargin{ArgCount};
      else
         error([mfilename '.m--Line colour specification defined more than once.']);
      end % if      
      
   end % if

end % for

% Use defaults for c and CMap if not specified.
if exist('CMap')==0,
   CMap = get(gcf,'colormap');
end % if

if exist('c')==0,
   if exist('z')==1,
      c = z;
   else
      c = [1:length(x)];
   end % if
end % if

% Index the supplied colourmap.
ColourIndex = round(interp1([min(c) max(c)],[1 length(CMap)],c));
keyboard
% Only make as many lines as there are separate colour segments.
GroupedIndex = groupneighbours(ColourIndex);
unique(ColourIndex)


NumLines = length(x)-1;
LineHndls = NaN * ones(NumLines,1);

for LineCount = 1:NumLines,

   if exist('z')==1,
      LineHndls(LineCount) = line([x(LineCount) x(LineCount+1)],[y(LineCount) y(LineCount+1)],[z(LineCount) z(LineCount+1)]);
   else
      LineHndls(LineCount) = line([x(LineCount) x(LineCount+1)],[y(LineCount) y(LineCount+1)]);
   end % if

   set(LineHndls(LineCount),'color',CMap(ColourIndex(LineCount),:));

end % for


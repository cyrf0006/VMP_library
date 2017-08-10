function groupMembers = groupneighbours(IndexIn)
%
% groupneighbours.m--Takes a vector of numbers (usually an index to another
% vector) and groups them into contiguous segments.
%
% Groupneighbours.m returns groupMembers, a cell array, each cell of which
% contains a group of next-door neighbours. A vector element with no other
% element adjacent to it will appear in its own cell.  
%
% Syntax: groupMembers = groupneighbours(IndexIn);
%
% e.g., groupMembers = groupneighbours([20 30 40 41 42 43 50 60 61 70 71 72 80])

% Kevin Bartlett (bartlettk@dfo-mpo.gc.ca) 11/1999
%------------------------------------------------------------------------------
% Test for development: groupMembers = groupneighbours([20 30 40 41 42 43 50 52 55 57 60 61 70 71 72 80 81])

% Modification, 2003-03-14: program very slow on large vectors. Speed up by
% looping over groups of neighbours, rather than individual points.

if isempty(IndexIn),
   groupMembers = [];
   return;
end % if

% Convert to column vector if necessary.
indexSize = size(IndexIn);

if ~any(indexSize==1)
   error([mfilename '.m--IndexIn must be a vector.'])
end % if

if indexSize(1)>indexSize(2)
   IndexIn = IndexIn';
   wasFlipped = 1;
else
   wasFlipped = 0;
end % if

% Sort the vector of integers.
IndexIn = sort(IndexIn);

% Test that there are no repeated values in the input vector.
if length(unique(IndexIn)) ~= length(IndexIn),
   error([mfilename '.m--IndexIn must not have repeated values.'])
end % if

% Find which points in the vector are NOT adjacent to each other. These
% points mark the divisions ("fences") between groups of neighbours.
FenceIndex = find(diff(IndexIn)>1);

% The start of each group follows each division.
GroupStartIndex = FenceIndex + 1;

% ...The first element in the input vector is always the start of its group.
GroupStartIndex = [1 GroupStartIndex];

% The end of each group is given by the fence index, with the addition of
% the last point in the input vector, which is always the end of its group.
GroupEndIndex = [FenceIndex length(IndexIn)];

% Initialise output variable.
NumGroups = length(GroupStartIndex);
groupMembers = cell(1,NumGroups);

% Find the members of each group.
for GroupCount = 1:NumGroups
   groupMembers{GroupCount} = IndexIn([GroupStartIndex(GroupCount):GroupEndIndex(GroupCount)]);
end % for

% Flip the index to have the same orientation as the input index.
if wasFlipped == 1
   groupMembers = groupMembers';
end % if

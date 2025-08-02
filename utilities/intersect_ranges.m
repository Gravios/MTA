function intersection = intersect_ranges(ranges1, ranges2)

numRanges1 = size(ranges1,1);
numRanges2 = size(ranges2,1);
v1 = zeros([numRanges1,1]);
v2 = zeros([numRanges2,1]);

intersection = [];

for ind1 = 1:numRanges1
    for ind2 = 1:numRanges2
        low = max( ranges1(ind1, 1), ranges2(ind2, 1) );
        high = min( ranges1(ind1, 2), ranges2(ind2, 2) );
        if low <= high %& ( ~isempty(low) & ~isempty(high) )
            intersection = [intersection;  low, high];
            v1(ind1) = 1;
            v2(ind2) = 1;
        end
    end
end
intersection = [];
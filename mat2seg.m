function [seg,segments] = mat2seg(mat,minCoverage,minSegment)
% TODO: maybe return segments as well
% Coverage is a scalar controlling how much % of the image we want to cover
if nargin < 2, minCoverage = 1; end
if nargin < 3, minSegment  = 0; end
assert(isscalar(minCoverage) && minCoverage > 0 && minCoverage <= 1, ...
    'minCoverage must be a scalar in (0,1]')
assert(isscalar(minSegment), 'minSegment must be scalar')

% Using this function assumes you have already grouped the medial points
% into branches. A "refined" MAT (using function refineMAT()) is not
% necessary, although it might lead to better results.
if ~isfield(mat,'branches') || isempty(mat.branches)
    mat.branches = groupMedialPoints(mat);
end

% Compute the depth contribution of each branch separately.
[H,W] = size(mat.depth);
numBranches = max(mat.branches(:));
depthBranch = zeros(H,W,numBranches);
for i=1:numBranches
    depthBranch(:,:,i) = mat2mask(mat.radius .* double(mat.branches == i),mat.scales);
end

% Segments are the areas covered by individual branches.
segments = double(depthBranch > 0);

% Sort by segment "importance", which is proportional to the area covered.
% Because of potential grouping errors, significant areas of the image may
% be covered by multiple segments, so we must take into account only the
% *new* pixels covered by each segment, by using this hack:
[~, idxSorted] = sort(sum(sum(segments)), 'descend');
segments = segments(:,:,idxSorted);
sumSeg   = cumsum(segments,3);
segments = ((segments - sumSeg) == 0) & (sumSeg > 0);
[areaSorted, idxSorted] = sort(sum(sum(segments)), 'descend');
segments = segments(:,:,idxSorted);

% Assign a different label to each segment. After sorting, the smaller the
% label, the larger the respective segment.
segments = bsxfun(@times, segments, reshape(1:numBranches,1,1,[]));

% Discard small segments
if minSegment > 0
    if minSegment < 1   % ratio of the min segment area over image area
        small = areaSorted/(H*W) < minSegment;
    elseif minSegment < H*W  % #pixels of min segment
        small = areaSorted < minSegment;
    else
        error('minSegment is larger than the size of the image')
    end 
    % If no segment satisfies the contraint, just use the largest segment
    if numel(small) == numel(areaSorted)
        small(1) = false;
    end
    segments(:,:,small) = [];
    areaSorted(small)   = [];
end

% Keep segments that cover at least (minCoverage*100) % of the image area.
if minCoverage < 1
    cumAreaSorted = cumsum(areaSorted)/(H*W);
    numSegmentsKeep = find(cumAreaSorted >= minCoverage, 1);
    if isempty(numSegmentsKeep)
        numSegmentsKeep = numel(cumAreaSorted);
        warning('%.1f%% coverage achieved (<%.1f%%)',...
            cumAreaSorted(numSegmentsKeep)*100,minCoverage*100)
    end
    segments = segments(:,:,1:numSegmentsKeep);
end
seg = max(segments,[],3);

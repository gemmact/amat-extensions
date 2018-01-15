function e = mat2edges(mat,r,s)
% Using this function assumes you have already grouped the medial points
% into branches. A "refined" MAT (using function refineMAT()) is not
% necessary, although it might lead to better results.
if ~isfield(mat,'branches')
    mat.branches = groupMedialPoints(mat);
end

% Default parameters for nonmaximum suppression.
if nargin < 2, r = 5; end
if nargin < 3, s = 5; end

branches = mat.branches;
radius   = mat.radius;
% Compute the depth contribution of each branch separately.
[H,W] = size(mat.depth);
numBranches = max(branches(:));
depthBranch = zeros(H,W,numBranches);
for i=1:numBranches
    depthBranch(:,:,i) = mat2mask(radius .* double(branches == i));
end

branchImportance = squeeze(sum(sum(depthBranch>0)));
[~,idxSorted] = sort(branchImportance,'descend');

% The edges are the points in the image that are covered by few disks
% belonging in the current branch AND by few disks belonging to OTHER
% branches. We use the edge strength of all the other branches to weigh the
% strength of the edges of the current branch.
d  = depthBranch;
dc = bsxfun(@minus,mat.depth,d); % complement of the depth for a branch
d  = 1-bsxfun(@rdivide,d, max(eps,max(max(d))));
dc = 1-bsxfun(@rdivide,dc,max(eps,max(max(dc))));
dp = d .* dc;

% Use the most important branches to compute the boundaries
topBranches = max(min(10,numBranches), ceil(numBranches/4));
e = prod(dp(:,:,idxSorted(1:topBranches)),3);

% Nonmaximum suppression (optional)
if nargin > 1
    e = edgeNms(single(e),r,s);
end


% Taken from Piotr Dollar's older version of structured edge detection code
% -------------------------------------------------------------------------
function E = edgeNms( E, r, s )
% -------------------------------------------------------------------------
% Suppress locations where edge is stronger in orthogonal direction
O = edgeOrient(E,r);
E1=padarray(E,[r+1 r+1],'replicate'); Dx=cos(O); Dy=sin(O);
[ht,wd]=size(E1); [cs,rs]=meshgrid(r+2:wd-r-1,r+2:ht-r-1);
for i=-r:r, if(i==0), continue; end
  cs0=i*Dx+cs; dcs=cs0-floor(cs0); cs0=floor(cs0);
  rs0=i*Dy+rs; drs=rs0-floor(rs0); rs0=floor(rs0);
  E2 = (1-dcs).*(1-drs) .* E1(rs0+0+(cs0-1)*ht);
  E2 = E2 + dcs.*(1-drs) .* E1(rs0+0+(cs0-0)*ht);
  E2 = E2 + (1-dcs).*drs .* E1(rs0+1+(cs0-1)*ht);
  E2 = E2 + dcs.*drs .* E1(rs0+1+(cs0-0)*ht);
  E(E*1.01<E2) = 0;
end
% suppress noisy estimates near boundaries
for r=1:s, E([r end-r+1],:,:)=E([r end-r+1],:,:)*(r-1)/s; end
for r=1:s, E(:,[r end-r+1],:)=E(:,[r end-r+1],:)*(r-1)/s; end

% -------------------------------------------------------------------------
function O = edgeOrient( E, r )
% -------------------------------------------------------------------------
% Compute very approximate orientation map from edge map. We have replaced
% convTri with the slightly slower equivalent code to reduce dependencies.
% E2=convTri(E,r); 
% ----------------------------------------------
if(r<=1), p=12/r/(r+2)-2; f=[1 p 1]/(2+p); r=1;
else f=[1:r r+1 r:-1:1]/(r+1)^2; end
E2 = padarray(E,[r r],'symmetric','both');
E2 = convn(convn(E2,f,'valid'),f','valid');
% ----------------------------------------------
f=[-1 2 -1];
Dx=conv2(E2,f,'same'); Dy=conv2(E2,f','same');
F=conv2(E2,[1 0 -1; 0 0 0; -1 0 1],'same')>0;
Dy(F)=-Dy(F); O=mod(atan2(Dy,Dx),pi);

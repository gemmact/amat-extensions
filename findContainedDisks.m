function [xc,yc,rc] = findContainedDisks(p,R,sz)
% Find center coordinates and respective radii of all disks that are FULLY
% contained in any given disk of radius R. Note that this function returns
% RELATIVE coordinates with respect to the disk's center.

% Center coordinates
xp = c(1); yp = c(2);
% Grid for a disk of radius R
[x,y] = meshgrid(-R:R,-R:R);
% Array centers of contained disks at each scale
dc = bsxfun(@le,x.^2 + y.^2,reshape(((R-1):-1:0).^2, 1,1,[]));
[yc,xc,rc] = ind2sub([2*R+1,2*R+1,R],find(dc));
yc = yc - R - 1 + yp; xc = xc - R - 1 + xp;
% If image size (sz) is given, remove invalid points that cross the boundary
if nargin == 3
    outOfLimits = xc < 1 | yc < 1 | xc > W | yc > H;
    xc(outOfLimits) = []; yc(outOfLimits) = []; rc(outOfLimits) = [];
end
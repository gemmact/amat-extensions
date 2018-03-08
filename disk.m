function d = disk(r,method)
% DISK Returnns a logical mask of size (2*r+1)x(2*r+1), in the shape of a
%   disk of radius r.
% 
%   d = DISK(r,method) controls the method used to construct the disk. The
%   method parameter can take two values:
% 
%   'simple' {default}
%   'alternative'
% 
%   Please check the source code for more information.
%
%   See also: meshgrid
% 
% Stavros Tsogkas <tsogkas@cs.toronto.edu>
% Last update: November 2016

if nargin < 2, method = 'simple'; end

r = double(r); % make sure r can take negative values
switch method
    case 'simple'
        [x,y] = meshgrid(-r:r, -r:r);
        d = x.^2 + y.^2 <= r^2;
    case {'alternative','alt'}
        r = r+1;
        [x,y] = meshgrid(-r:r, -r:r);
        d = x.^2 + y.^2 < r^2;
        d = d(2:(end-1),2:(end-1));
    case 'quadrant1'
        [x,y] = meshgrid(-r:r, -r:r);
        d = x.^2 + y.^2 <= r^2;
        d = d(1:r+1,1:r+1);
    case 'quadrant2'
        [x,y] = meshgrid(-r:r, -r:r);
        d = x.^2 + y.^2 <= r^2;
        d = d(1:r+1,r+1:end);
    case 'quadrant3'
        [x,y] = meshgrid(-r:r, -r:r);
        d = x.^2 + y.^2 <= r^2;
        d = d(r+1:end,1:r+1);
    case 'quadrant4'
        [x,y] = meshgrid(-r:r, -r:r);
        d = x.^2 + y.^2 <= r^2;
        d = d(r+1:end,r+1:end);
    otherwise
        error('Method not supported')
end

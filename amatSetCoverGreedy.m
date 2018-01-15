%% Setup global parameters and preprocess image
R = 50; % #scales
B = 32; % #bins
errorType = 'se';

% imgRGB = rgb2gray(im2double(imresize(imread('google.jpg'), [128 128], 'nearest'))); 
% imgRGB = rgb2gray(im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/66075.jpg'), [128 128], 'nearest'))); 
% imgRGB = im2double(imresize(imread('/home/tsogkas/datasets/BSDS500/images/train/35070.jpg'), [128 128], 'nearest')); 
imgRGB = rgb2gray(im2double(imresize(imread('/home/tsogkas/datasets/AbstractScenes_v1.1/RenderedScenes/Scene0_9.png'),0.5)));
[H,W,C] = size(imgRGB);
imgLab = rgb2labNormalized(imgRGB);
% imgClustered = clusterImageValues(imgLab, 5); % simplify input

%% Construct filters, calculate perimeneters and disk areas
filters = cell(R,1); for r=1:R, filters{r} = disk(r); end

%% Compute encodings f(D_I(x,y,r)) at every point.
% mrgb = imageEncoding(imgRGB,filters,'average'); % used for the reconstruction error
mlab = imageEncoding(imgLab,filters,'average'); % used for the reconstruction error
% hlab = imageEncoding(binImage(imgLab,B),filters,'hist',B); % used for the maximality error;
% htex = imageEncoding(textonMap(imgRGB,B),filters,'hist',B);
htex = []; % NO TEXTURE FOR NOW

%% Compute reconstruction and maximality errors at all points and scales
%  Compute reconstruction error by using the dssim (structural
%  dissimilarity metric in the RGB color space, or a (nm)se metric on the
%  perceptually linear Lab color space. 
if strcmp(errorType,'dssim')
    reconstructionError = imageError(imgRGB,mrgb,filters,'dssim');
else
    luminanceError = imageError(imgLab(:,:,1), mlab(:,:,1,:),filters, errorType);
    if ismatrix(imgLab)
        reconstructionError = luminanceError;
    else
        colorError = imageError(imgLab(:,:,2), mlab(:,:,2,:), filters, errorType) + ...
                     imageError(imgLab(:,:,3), mlab(:,:,3,:), filters, errorType);
        reconstructionError = (luminanceError + colorError)/2;
    end
end
%  The maximality error that encourages the selection of maximal disks
maximalityError = zeros(H,W,R);
h = cat(3,hlab,htex);
for r=1:R-1
    dr = ceil(r/(2+sqrt(6))); % dA >= 0.5A(r)
%     dr = ceil(r/(1+sqrt(2))); % dA == A(r)
%     dr = 1;
    h1 = h(:,:,:,:,r);
    h2 = h(:,:,:,:,min(r+dr,R))-h(:,:,:,:,r);
    h1 = bsxfun(@rdivide,h1,sum(h1,4));
    h2 = bsxfun(@rdivide,h2,sum(h2,4));
    maximalityError(:,:,r) = sum(...
        histogramDistance(h2(:,:,1:end-1,:), h1(:,:,1:end-1,:),'chi2-gaussian',0.2),3) + ...
        2*histogramDistance(h2(:,:,end,:), h1(:,:,end,:),'chi2');        
end
maximalityError = max(0,1-maximalityError);
% Fix boundary conditions for both  terms (image boundaries crosses)
for r=1:R
    reconstructionError(1:r,:,r)         = inf;
    reconstructionError(:,1:r,r)         = inf;
    reconstructionError(end-r+1:end,:,r) = inf;
    reconstructionError(:,end-r+1:end,r) = inf;
    maximalityError(1:r,:,r)             = 1;
    maximalityError(:,1:r,r)             = 1;
    maximalityError(end-r+1:end,:,r)     = 1;
    maximalityError(:,end-r+1:end,r)     = 1;
end
maximalityError(:,:,1:2) = 1;

%% Compute fully contained disks for a point in the image
% [xc,yc,rc] = findContainedDisks(r);
% This is the grid used for all scales from 1 to R. Centers for contained
% disks at all scales can be deduced by using this grid and adding the
% necessary offset depending on the position of the current point.
[x,y] = meshgrid(-R:R,-R:R);
% Centers of all contained disks of scales 1,...,R.
dc = bsxfun(@le,x.^2 + y.^2,reshape(((R-1):-1:0).^2, 1,1,[]));
% We want to compute the sum of square differences of the encoding (mean) 
% of each point in the image from the encodings of all the included disks. 
% EDIT: WE CANNOT USE THE EXACT SAME SIMPLIFIED FORMULA AS IN imageError() because
% that derivation assumed that g = mean(I) for some patch, which is not the
% case here, as each of the I values is a mean value of the original image
% itself. We can still use convolutions though.
% Precompute areas. A(:,:,end-r+1) is the map with the pixels covered by a
% disk of radius r for the whole image domain.
A = ones(H,W,R);
for r=1:R
    D = dc(:,:,r); D = cropImageBox(D,mask2bbox(D)); 
    A(:,:,r) = conv2(A(:,:,r),double(D),'same');
end
A = cumsum(A,3,'reverse');
reconstructionError = zeros(H,W,R);
enc  = mlab; enc2 = enc.^2;
for r=1:R % compute consensus for all scales
    dcsubset = dc(:,:,end-r+1:end);
    % for a given disk at scale r, consider all radii of contained disks
    for i=1:size(dcsubset,3) 
        D = dcsubset(:,:,i); D = double(cropImageBox(D,mask2bbox(D))); 
        reconstructionError(:,:,r) = ... 
            reconstructionError(:,:,r) + ...
            conv2(enc2(:,:,1,i),D,'same') - ...
            conv2(enc(:,:,1,i), D,'same') .* 2 .* enc(:,:,1,r);
    end
%     consensusScores(:,:,r) = consensusScores(:,:,r) + A * enc2(:,:,1,r);
    reconstructionError(:,:,r) = reconstructionError(:,:,r) + A(:,:,end-r+1) .* enc2(:,:,1,r);
end
reconstructionError = max(0,reconstructionError);

% Compute maximality scores using the mean value consensus
% Create ring masks of outer rings
maximalityError = zeros(H,W,R);
for r=1:R
    dr = ceil(r/(2+sqrt(6))); % dA >= 0.5A(r)
    mask = double(ring(r,r+dr));
    maximalityError(:,:,r) = 1-(conv2(imgRGB.^2,mask,'same') + ...
        nnz(mask)*enc2(:,:,1,r) - 2 .* enc(:,:,1,r) .* conv2(imgRGB,mask,'same'))/nnz(mask);
end 
maximalityError = max(0,maximalityError);

% Fix boundary conditions for both  terms (image boundaries crosses)
for r=1:R
    reconstructionError(1:r,:,r)         = inf;
    reconstructionError(:,1:r,r)         = inf;
    reconstructionError(end-r+1:end,:,r) = inf;
    reconstructionError(:,end-r+1:end,r) = inf;
    maximalityError(1:r,:,r)             = inf;
    maximalityError(:,1:r,r)             = inf;
    maximalityError(end-r+1:end,:,r)     = inf;
    maximalityError(:,end-r+1:end,r)     = inf;
end


% We can obtain the centers of contained disks at all scales for a disk of
% radius r, by using dc(:,:,end-r+1:end)
% [yc,xc,rc] = ind2sub([2*R+1,2*R+1,R],find(dc));
% contains = false(H,W,R); contains(sub2ind([H,W,R],yc,xc,rc)) = true;
% Adjust indices (which are always positive) to local grid coordinates and
% add offset to center at the current point.
% [x,y] = meshgrid(1:W,1:H);
% xp = 100; yp = 100;
% yc = yc - R - 1 + yp; xc = xc - R - 1 + xp;
% outOfLimits = xc < 1 | yc < 1 | xc > W | yc > H;
% xc(outOfLimits) = []; yc(outOfLimits) = []; rc(outOfLimits) = [];



%% Greedy approximation of the weighted set cover problem associated with AMAT
% Initializations
amat.input          = reshape(imgLab, H*W, C);
amat.reconstruction = reshape(amat.input, H*W, C);
amat.axis           = zeros(H,W,C);
amat.radius         = zeros(H,W);
amat.depth          = zeros(H,W); % #disks points(x,y) is covered by
amat.covered        = false(H,W);
amat.price          = zeros(H,W); % error contributed by each point

%% Error balancing and visualization of top (low-cost) disks
% -------------------------------------------------------------------------
% Disk cost: refers to the cost contributed to the total absolute cost of 
% placing the disk in the image.
% -------------------------------------------------------------------------
% Disk cost effective: refers to the "normalized" cost of the disk, divided
% by the number of *new* pixels it covers when it is added in the solution.
% -------------------------------------------------------------------------
% numNewPixelsCovered: number of NEW pixels that would be covered by each
% disk if it were to be added in the solution. This number can be easily 
% computed using convolutions, while taking into account the fact that 
% some disks exceed image boundaries.
% -------------------------------------------------------------------------
numNewPixelsCovered = ones(H,W,R);
for r=1:R
    numNewPixelsCovered(:,:,r) = conv2(numNewPixelsCovered(:,:,r), ...
        double(filters{r}),'same');
end
% diskCostEffective = min(1,sqrt(reconstructionError ./ numNewPixelsCovered) + bsxfun(@plus,maximalityError,reshape(T./(1:R),1,1,[])));
diskCostEffective = reconstructionError ./ numNewPixelsCovered + maximalityError;
% diskCostEffective = reconstructionError + maximalityError;
% Sort costs in ascending order and visualize top disks.
[sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
top = 1e2;
[yy,xx,rr] = ind2sub([H,W,R], indSorted(1:top));
figure(1); imshow(reshape(amat.input,H,W,[])); 
viscircles([xx,yy],rr,'Color','k','LineWidth',0.5);
viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
title(sprintf('W: Top-%d disks, B: Top-1 disk',top))
errorBackup = reconstructionError;


%% Run the greedy algorithm
[x,y]   = meshgrid(1:W,1:H);
% [xr,yr] = meshgrid(1:R,1:R);
% containedDisks = bsxfun(@le,xr.^2 + yr.^2,reshape(((R-1):-1:0).^2, 1,1,[]));
% coveringDisks  = bsxfun(@le,xr.^2 + yr.^2,reshape((1:R).^2, 1,1,[]));
% [ycov,xcov,rcov] = ind2sub([2*R+1,2*R+1,R],find(dc));
% ycov = ycov-R-1; xcov = xcov-R-1;
f = mlab;
while ~all(amat.covered(:))
    % Find the most cost-effective set at the current iteration
    [minCost, indMin] = min(diskCostEffective(:));
    % D is the set of points covered by the selected disk
    [yc,xc,rc] = ind2sub([H,W,R], indMin);
    distFromCenterSquared = (x-xc).^2 + (y-yc).^2;
    D = distFromCenterSquared <= rc^2;
    newPixelsCovered = D & ~amat.covered;
    assert(~all(newPixelsCovered(:)))
    
    % Update AMAT
    reconstructedDisk = repmat(reshape(f(yc,xc,:,rc),[1 C]), [nnz(newPixelsCovered),1]);
    amat.price(newPixelsCovered) = sum(( ...
        amat.reconstruction(newPixelsCovered,:) - reconstructedDisk ).^2,2);
    amat.covered(D) = true;
    amat.depth(D) = amat.depth(D) + 1;
    amat.reconstruction(newPixelsCovered,:) = reconstructedDisk;
    amat.axis(yc,xc,:) = f(yc,xc,:,rc);
    amat.radius(yc,xc) = rc;

    % Find how many of the newPixelsCovered are covered by other disks in
    % the image and subtract the respective counts from those disks.
    [yy,xx] = find(newPixelsCovered);
    xmin = min(xx); xmax = max(xx);
    ymin = min(yy); ymax = max(yy);
    newPixelsCovered = double(newPixelsCovered);
    priceMap = amat.price .* newPixelsCovered;
    amat.reconstruction = reshape(amat.reconstruction,H,W,[]); % reshape for convolutions
    for r=1:R
        xxmin = max(xmin-r,1); yymin = max(ymin-r,1);
        xxmax = min(xmax+r,W); yymax = min(ymax+r,H);
        numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) = ...
            numNewPixelsCovered(yymin:yymax,xxmin:xxmax, r) - ...
            conv2(newPixelsCovered(yymin:yymax,xxmin:xxmax),double(filters{r}),'same');
        reconstructionError(yymin:yymax,xxmin:xxmax, r) = ...
            reconstructionError(yymin:yymax,xxmin:xxmax, r) - ...
            conv2(priceMap(yymin:yymax,xxmin:xxmax),double(filters{r}),'same');
        % TODO: use conv2 with 'valid' parameter to contain and speed up?
        % We must do this "manually" if we want to use imageEncoding() by
        % using the larger window [xxmin,yymin,xxmax,yymax] to compute the
        % convolutions and then crop the [xmin,ymin,xmax,ymax] part.
        tmp = imageEncoding(amat.reconstruction(yymin:yymax,xxmin:xxmax,:),filters(r));
        f(ymin:ymax,xmin:xmax, :, r) = ...
            tmp((ymin-yymin+1):(ymax-yymin+1),(xmin-xxmin+1):(xmax-xxmin+1), :);
    end
    amat.reconstruction = reshape(amat.reconstruction,H*W,[]); % reshape back

    
    % We do this to correct precision errors introduced by subtracting the
    % pixel prices from the reconstruction error (these would normally be
    % zero but sometimes they take very small negative values).
    reconstructionError = max(0,reconstructionError); 
    
    % Update errors. NOTE: the diskCost for disks that have
    % been completely covered (e.g. the currently selected disk) will be
    % set to inf or nan, because of the division with numNewPixelsCovered
    % which will be zero (0) for those disks. 
    diskCostEffective = reconstructionError ./ numNewPixelsCovered + maximalityError;
    diskCostEffective(isnan(diskCostEffective)) = inf;
    assert(allvec(numNewPixelsCovered(yc,xc, 1:rc)==0))

    
    % Visualize progress
    if 1
        % Sort costs in ascending order to visualize updated top disks.
        [sortedCosts, indSorted] = sort(diskCostEffective(:),'ascend');
        [yy,xx,rr] = ind2sub([H,W,R], indSorted(1:top));
        subplot(221); imshow(reshape(amat.input, H,W,[])); 
        viscircles([xc,yc],rc, 'Color','k','EnhanceVisibility',false);
        title('Selected disk');
        subplot(222); imshow(bsxfun(@times, reshape(amat.input,H,W,[]), double(~amat.covered))); 
        viscircles([xx,yy],rr,'Color','w','EnhanceVisibility',false,'Linewidth',0.5); 
        viscircles([xx(1),yy(1)],rr(1),'Color','b','EnhanceVisibility',false); 
        viscircles([xc,yc],rc,'Color','y','EnhanceVisibility',false); 
        title(sprintf('K: covered %d/%d, W: Top-%d disks,\nB: Top-1 disk, Y: previous disk',nnz(amat.covered),H*W,top))
        subplot(223); imshow(amat.axis); title('A-MAT axes')
        subplot(224); imshow(amat.radius,[]); title('A-MAT radii')
        drawnow;
    end
    disp(nnz(~amat.covered))
end
amat.reconstruction = reshape(amat.reconstruction,H,W,C);
amat.input = reshape(amat.input,H,W,C);

%% Visualize results
figure(3); clf;
subplot(221); imshow(amat.axis); title('Medial axes');
subplot(222); imshow(amat.radius,[]); title('Radii');
subplot(223); imshow(amat.input); title('Original image');
subplot(224); imshow(amat.reconstruction); title('Reconstructed image');

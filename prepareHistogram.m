function ii = prepareHistogram(im, B)
%Calculates integral images for each bin

[H, W, C] = size(im);

w = (1/B);
ii = zeros(H, W, B);
imgBinned = binImage(im,B);
% ii(:, :, 1) = cumsum(cumsum(im == 0), 2);
for j = 1:B
%     ij = (im > (j-1)*w & im <= j*w);
    ii(:, :, j) = cumsum(cumsum( imgBinned == j), 2);
    
end
end
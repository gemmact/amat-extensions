function ii = prepareHistogram(im, B)
%Calculates integral images for each bin

[H, W, C] = size(im);

w = ceil(255/B);
ii = zeros(H, W, B);
ii(:, :, 1) = cumsum(cumsum(im == 0), 2);
for j = 1:B
    ij = (im > (j-1)*w & im <= j*w);
    ii(:, :, j) = cumsum(cumsum(ij), 2);
    
end
end
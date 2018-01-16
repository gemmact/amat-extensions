function h = createHistogram(ii, x1, y1, P)
    %use integral images to create normalized histogram
    [H, W, B] = size(ii);

    x2 = x1 + P-1;
    y2 = y1 + P-1;
    
    
    h = zeros(1, B);
    for i = 1:B
        h(i) = ii(x2, y2, i) + ii(x1, y1, i) - ii(x1, y2, i) - ii(x2, y1, i);
    end
    
    %h = h/(P*P);

end
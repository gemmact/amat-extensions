function h = createHistogram(ii, x1, y1, P)
    %use integral images to create normalized histogram
    [H, W, B] = size(ii);

    x2 = x1 + P-1;
    y2 = y1 + P-1;
    
    
    h = zeros(1, B);
    for i = 1:B
        h(i) = ii(y2, x2, i) + ii(y1, x1, i) - ii(y1, x2, i) - ii(y2, x1, i);
    end
    
    %h = h/(P*P);

end
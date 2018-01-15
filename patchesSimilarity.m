
function h = patchesSimilarity(im)
global B

patchSize = 32;
B = 32;

global h1

im = imread(im);
im = rgb2labNormalized(im);


[H,W,C] = size(im);
RI = imref2d(size(im));
RI.XWorldLimits = [0 H];
RI.YWorldLimits = [0 W];
fh = figure(1); imshow(im, RI, 'InitialMagnification', 'fit');
set(fh, 'WindowButtonDownFcn', @firstPoint);
set(fh, 'WindowButtonMotionFcn', @changePoint);

global ii1 ii2 ii3 ii4

ii1 = prepareHistogram(im(:, :, 1), B);
ii2 = prepareHistogram(im(:, :, 2), B);
ii3 = prepareHistogram(im(:, :, 3), B);

img = rgb2gray(im2double(im)); 
texture = textonMap(img, B);
ii4 = prepareHistogram(texture, B);



end

function firstPoint(fh,~)
        global h1 B ii1 ii2 ii3 ii4
        im = imread('9170.jpg');
        [H,W,C] = size(im);
        
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        
        RI = imref2d(size(im));
        RI.XWorldLimits = [0 H];
        RI.YWorldLimits = [0 W];
        fh = figure(1); imshow(im, RI, 'InitialMagnification', 'fit');
        
        
        P = 32; %patchsize
        x1 = x - floor(P/2);
        x2 = x + ceil(P/2);
        y1 = y - floor(P/2);
        y2 = y + ceil(P/2);
        
        if (x1>= 0 && x2<= H && y1>=0 && y2<=W)
            
            rectangle('Position', [x1, y1, P, P], 'EdgeColor', 'r');

            h1 = zeros(84, B);
            for k = 2.^(0:2)
                patchSize = P/k;
                for i = x1:x2-patchSize+1
                    for j = y1:y2-patchSize+1
                        h1(4*k-3, :) = createHistogram(ii1, i, j, patchSize);
                        h1(4*k-2, :) = createHistogram(ii2, i, j, patchSize);
                        h1(4*k-1, :) = createHistogram(ii3, i, j, patchSize);
                        h1(4*k, :) = createHistogram(ii4, i, j, patchSize);
                        
                    end
                end
            end
        end
            
                   
            %reference_patch = im((x1-patchSize/2 : x1 + patchSize/2), (y1-patchSize/2 : y1 + patchSize/2), :)
            %im((x1-patchSize/2 : x1 + patchSize/2), (y1-patchSize/2 : y1 + patchSize/2), :) = zeros(patchSize+1, patchSize+1, 3);
            %figure(2)

            %B = 5;
            %f1 = patchEncoding(binImage(reference_patch, B), 'hist-normalized', B); 
            
        
end



function changePoint(fh,~)
        global h1 B ii1 ii2 ii3 ii4
        im = imread('9170.jpg');
        [H,W,C] = size(im);
        P = 32;
        
        
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        
        x1 = x - floor(P/2);
        x2 = x + ceil(P/2);
        y1 = y - floor(P/2);
        y2 = y + ceil(P/2);
        
        P = 32; %patchsize
        
        
        if (x1>= 0 && x2<= H && y1>=0 && y2<=W)
            
            %r = rectangle('Position', [x1, y1, P, P]);
            
            h2 = zeros(84, B);
            for k = 2.^(0:2)
                patchSize = P/k;
                for i = x1:x2-patchSize+1
                    for j = y1:y2-patchSize+1
                        h2(4*k-3, :) = createHistogram(ii1, i, j, patchSize);
                        h2(4*k-2, :) = createHistogram(ii2, i, j, patchSize);
                        h2(4*k-1, :) = createHistogram(ii3, i, j, patchSize);
                        h2(4*k, :) = createHistogram(ii4, i, j, patchSize);
                    end
                end
            end



                %second_patch = im((x - patchSize/2 : x + patchSize/2), (y-patchSize/2 : y + patchSize/2), :);
                %%imshow(reference_patch);

                %f1 = patchEncoding(binImage(reference_patch, B), 'hist', B) 
                %B = 5;
                %f2 = patchEncoding(binImage(second_patch, B), 'hist-normalized', B);
                
                


                dist = histogramDistance(h1/norm(h1), h2/norm(h2), 'chi2');
                
                
                
        end
end


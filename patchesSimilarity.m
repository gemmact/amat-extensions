
function h = patchesSimilarity(im)
global B imRGB

imRGB = imread(im);

patchSize = 32;
B = 32;

global x1ref y1ref


imLab = rgb2labNormalized(imRGB);

x1ref = -32;
y1ref = -32;



[H,W,C] = size(imRGB);
RI = imref2d(size(imRGB));
RI.XWorldLimits = [1 W];
RI.YWorldLimits = [1 H];
fh = figure(1); 
subplot(3,6,[1 15])
imshow(imRGB, RI, 'InitialMagnification', 'fit');
set(fh, 'WindowButtonDownFcn', @firstPoint);
set(fh, 'WindowButtonMotionFcn', @changePoint);

global ii1 ii2 ii3 ii4

ii1 = prepareHistogram(imLab(:, :, 1), B);
ii2 = prepareHistogram(imLab(:, :, 2), B);
ii3 = prepareHistogram(imLab(:, :, 3), B);

img = rgb2gray(im2double(imRGB)); 
texture = textonMap(img, B);
ii4 = prepareHistogram(texture, B);



end




function firstPoint(fh,~)
        global h1 B ii1 ii2 ii3 ii4 x1ref y1ref imRGB
        %im = imread(im);
        im = imRGB;
        [H,W,C] = size(im);
        
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        
        RI = imref2d(size(im));
        RI.XWorldLimits = [1 W];
        RI.YWorldLimits = [1 H];
         
        
        subplot(3,6,[1 15])
        imshow(im, RI, 'InitialMagnification', 'fit');
        
        
        
        P = 32; %patchsize
        %x = 50;
        %y = 50;
        x1 = x - floor(P/2);
        x2 = x + ceil(P/2) - 1;
        y1 = y - floor(P/2);
        y2 = y + ceil(P/2) - 1;
        
        x1ref = x1;
        y1ref = y1;
       
        
        if (x1>= 1 && x2<= W && y1>=1 && y2<=H)
            
            rectangle('Position', [x1, y1, P, P], 'EdgeColor', 'r');

            h1 = zeros(84, B);
            %h1 = zeros(4,B);
            l = 1;
            for k = 2.^(0:2)
                patchSize = P/k;
                for i = x1:patchSize:x2-patchSize+1
                    for j = y1:patchSize:y2-patchSize+1
                        
                        h1(4*l-3, :) = createHistogram(ii1, i, j, patchSize);
                        h1(4*l-2, :) = createHistogram(ii2, i, j, patchSize);
                        h1(4*l-1, :) = createHistogram(ii3, i, j, patchSize);
                        h1(4*l, :) = createHistogram(ii4, i, j, patchSize);
                        l = l + 1;
                        
                        
                    end
                end
            end
        end
            
                   
            reference_patch = im((y1 : y2), (x1 : x2), :);
            %figure(1)
            subplot(3,6,4);
            imshow(reference_patch);
            subplot(3,6,10);
            imshow(reference_patch);
            hold on;
            [X,Y]=meshgrid(0:16:33);
            plot(X,Y,'k');
            plot(Y,X,'k'); axis off
            
            
            subplot(3,6,16);
            imshow(reference_patch);
            hold on;
            [X,Y]=meshgrid(0:8:33);
            plot(X,Y,'k');
            plot(Y,X,'k'); axis off
            

            %B = 5;
            %f1 = patchEncoding(binImage(reference_patch, B), 'hist-normalized', B); 
            
        
end



function changePoint(fh,~)
        global h1 B ii1 ii2 ii3 ii4 x1ref y1ref imRGB
        im = imRGB;
        [H,W,C] = size(im);
        P = 32;
        
        RI = imref2d(size(im));
        RI.XWorldLimits = [1 W];
        RI.YWorldLimits = [1 H];
         
        subplot(3,6,[1 15])
        imshow(im, RI, 'InitialMagnification', 'fit');
        rectangle('Position', [x1ref, y1ref, P, P], 'EdgeColor', 'r');
        
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        
        x1 = x - floor(P/2);
        x2 = x + ceil(P/2) - 1;
        y1 = y - floor(P/2);
        y2 = y + ceil(P/2) - 1;
        
        P = 32; %patchsize
        
        
        if (x1>= 1 && x2<= W && y1>=1 && y2<=H)
            
            r = rectangle('Position', [x1, y1, P, P]);
            
            if (x1ref >= 1 && x1ref <= W-P+1 && y1ref >= 1 && y1ref <= H-P+1)
            
                h2 = zeros(84, B);
                %h2 = zeros(4,B);
                l = 1;
                for k = 2.^(0:2)
                    patchSize = P/k;
                    for i = x1:patchSize:x2-patchSize+1
                        for j = y1:patchSize:y2-patchSize+1

                            h2(4*l-3, :) = createHistogram(ii1, i, j, patchSize);
                            h2(4*l-2, :) = createHistogram(ii2, i, j, patchSize);
                            h2(4*l-1, :) = createHistogram(ii3, i, j, patchSize);
                            h2(4*l, :) = createHistogram(ii4, i, j, patchSize);
                            l = l + 1;

                        end
                    end
                end




                    second_patch = im((y1 : y2), (x1 : x2), :);
                    subplot(3,6,6)
                    imshow(second_patch);
                    subplot(3,6,12)
                    imshow(second_patch);
                    hold on;
                    [X,Y]=meshgrid(0:16:33);
                    plot(X,Y,'k');
                    plot(Y,X,'k'); axis off


                    subplot(3,6,18)
                    imshow(second_patch);
                    hold on;
                    [X,Y]=meshgrid(0:8:33);
                    plot(X,Y,'k');
                    plot(Y,X,'k'); axis off
                    %%imshow(reference_patch);

                    %f1 = patchEncoding(binImage(reference_patch, B), 'hist', B) 
                    %B = 5;
                    %f2 = patchEncoding(binImage(second_patch, B), 'hist-normalized', B);


                    dist = histogramDistance(h1/norm(h1), h2/norm(h2), 'chi2')

                    subplot(3,6,5, 'replace')
                    [X,Y]=meshgrid(0:1);
                    hold on;
                    plot(X,Y,'k');
                    plot(Y,X,'k'); axis off;
                    text(0, 0.5, num2str(1 - sum(dist(1:4))/4));


                    subplot(3,6,11, 'replace')
                    [X,Y]=meshgrid(0:2);
                    hold on;
                    plot(X,Y,'k');
                    plot(Y,X,'k'); axis off
                    l = 5;
                    for i = 0:1
                        for j = 1:-1:0
                            text(i, j+0.5, num2str(1 - sum(dist(l:l+3))/4));
                            l = l + 4;
                        end
                    end


                    subplot(3,6,17, 'replace')
                    [X,Y]=meshgrid(0:4);
                    hold on;
                    plot(X,Y,'k');
                    plot(Y,X,'k'); axis off
                    l = 21;
                    for i = 0:3
                        for j = 3:-1:0
                            text(i, j+0.5, num2str(1 - sum(dist(l:l+3))/4));
                            l = l + 4;
                        end
                    end


                    w1 = 0.125*ones(1,4);
                    w2 = 0.01875*ones(1,16);
                    w3 = 0.003125*ones(1,64);

                    %h = 1- sum(dist)/length(dist)

                    h = 1 - w1*dist(1:4) + w2*dist(5:20) + w3*dist(21:84)
            end
        end
end

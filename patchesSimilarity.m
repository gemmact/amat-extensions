function fh = patchesSimilarity(imgRGB)
% DRAWDISKONFIGUREINTERACTIVE Draw reconstructed disk patches on the 
%   original image and display useful information. This tool allows the 
%   user to interactively change various parameters, such as the size of
%   the disk, the type of the error function used, and the method used for
%   summarizing an image patch. The supported functions include:
% 
%   [hover mouse over figure]: change the coordinates of the disk center.
%   [left click]: change the error type used ({'se'},'mse','nmse','rse',
%                 'rmse','nrmse','dssim').
%   [right click]:  change the encoding method ({'average'},'hist').
%   [middle click]: change the number of bins used for the histogram
%                   computations.
%   [scroll wheel]: up/down (decrease/increase) the radius of the disk.
%
%   See also: patchEncoding
% 
%   Stavros Tsogkas <tsogkas@cs.toronto.edu>
%   Last update: November 2016

% Default parameters
r = 15;
numBins = 32; % used for histogram encodings
methods.error    = {'se','mse','nmse','rse','rmse','nrmse','dssim'};
methods.encoding = {'average','hist-smirnov','hist-expectation','hist-mode'};
encodingType  = 'hist-normalized';
errorType = 'se';
errorCounter = find(strcmp(methods.error, errorType));
encodingCounter = find(strcmp(methods.encoding, encodingType));

% Plot figure and set callbacks
fh = figure(1);
subplot(3,6,[1 15]);
imshow(imgRGB);
im = imgRGB;
set(fh, 'WindowButtonMotionFcn', @changePoint);
set(fh, 'WindowButtonDownFcn',   @changeReferencePatch);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
[H,W,C] = size(imgRGB);
imgGray = rgb2gray(imgRGB); 
imgTexture = textonMap(imgGray, numBins);
size(imgTexture)
%imgTexture    = reshape(imgTexture,H*W,[]);
%tmap    = textonMap(imgRGB, numBins); 
%tmap    = reshape(tmap,H*W,[]);
imgRGB  = reshape(imgRGB, [], C);
imgLab  = rgb2labNormalized(imgRGB);
refPatch= zeros(2*r+1);
refxy   = [];
encPatch = zeros((C+1)*(1 + r + 1 + (r + 1)/4), numBins);
encRefPatch = zeros((C+1)*(1 + r + 1 + (r + 1)/4), numBins);
isReferenceFrameSelected = false;
selectingNewReference = false;
[xx,yy] = meshgrid(1:W,1:H);
w1 = 0.05*ones(1,4);
w2 = 0.01875*ones(1,16);
w3 = 0.0078125*ones(1,64);
showHomogeneityHeatmap()


    function drawSquare(fh)
        % Get point coordinates and check for validity
        fh = figure(1);
        subplot(3,6,[1 15]);
        imshow(im);
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        if x < 1 || x > W || y < 1 || y > H
            title('You are outside of the figure'); drawnow; return
        end
        if x-r < 1 || x+r > W || y-r < 1 || y+r > H
            title('Square crosses the image boundary'); drawnow; return
        end
        
        % Square logical mask
        D = abs(xx-x) <= r & abs(yy-y) <= r;

        % The dssim metric should be used on RGB data
        if strcmp(errorType, 'dssim')
            imgPatch = imgRGB(D,:);
        else
            imgPatch = imgLab(D,:);
        end
        
        if isReferenceFrameSelected
            figure(1)
            % Keep plotting the reference frame
            rectangle('Position', [refxy(1)-r, refxy(2)-r, 2*r, 2*r], 'EdgeColor', 'g');
            
            % Encode patch and subpatches
            encPatch = subpatchesEncode(x, y);
            

            % Disable annoying docking error that clutters the command line
            if strcmp(fh.WindowStyle, 'docked')
                warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
            end
            
            
            % Display image and compare/save patches
            if selectingNewReference
                
                %save reference patch
                refPatch = reshape(imgPatch, [2*r + 1, 2*r + 1, 3]);  
                
                % Encode patch and subpatches
                encRefPatch = subpatchesEncode(x, y);
                
                %Create heatmap of similarity values between the reference
                %patch and the other ones
                
                showSimilarityHeatmap(encRefPatch);
                
                %plot patch of reference and subpatches to study
                figure(1)
                
                subplot(3,6,4);
                imshow(refPatch);
                
                subplot(3,6,10);
                displaySubpatches(refPatch, r + 1);
                
                subplot(3,6,16);
                displaySubpatches(refPatch, floor(r/2) + 1);
                
            else
                %save patch to compare with reference one
                comparePatch = reshape(imgPatch, [2*r + 1, 2*r + 1, 3]);
 
                %plot square we are comparing to reference one
                rectangle('Position', [x-r, y-r, 2*r, 2*r], 'EdgeColor', 'k');
                
                %plot patch and subpatches to study
                subplot(3,6,6)
                imshow(comparePatch);
                    
                subplot(3,6,12)
                displaySubpatches(comparePatch, r + 1);

                subplot(3,6,18)
                displaySubpatches(comparePatch, floor(r/2) + 1);
                
                %compare patches
                dist = histogramDistance(encRefPatch, encPatch, 'chi2');
                
                
                %plot similarity scores between histograms
                
                subplot(3,6,5, 'replace')
                displayScores(0, 1, dist);

                subplot(3,6,11, 'replace')
                displayScores(1, 5, dist);

                subplot(3,6,17, 'replace')
                displayScores(3, 21, dist);
                    
                    
                %Calculate weighted distance between the two patches        



                h = 1 - (w1*dist(1:4) + w2*dist(5:20) + w3*dist(21:84))
                
                

            end
            
            
            drawnow;
        else
            % Simply plot square frame of the patch
            rectangle('Position', [x-r, y-r, 2*r, 2*r], 'EdgeColor', 'k');
        end        
    end


    function changePoint(fh,~)
        drawSquare(fh);
    end

    function changeRadius(fh,callbackData)
        r = min(min(H,W)/2, max(1, r + callbackData.VerticalScrollCount));
        drawSquare(fh);
    end
    
    function changeNumBins(fh)
        validInput = false;
        dlgTitle = 'Change number of histogram bins';
        while ~validInput
            answer = inputdlg('Enter number of bins:',dlgTitle);
            if isempty(answer)
                validInput = true; % keep previous nBins
            else
                answer = answer{1};
                answer = str2double(answer);
                if isempty(answer) || answer <= 0
                    dlgTitle = 'Invalid input! #bins must be a positive scalar.';
                else
                    numBins = answer;
                    validInput = true;
                end
            end
        end
        drawSquare(fh);
    end

    function changeReferencePatch(fh,~)                
        if strcmp(fh.SelectionType, 'normal')
            x = round(fh.CurrentAxes.CurrentPoint(1,1));
            y = round(fh.CurrentAxes.CurrentPoint(1,2));
            refxy = [x,y];
            isReferenceFrameSelected = true;
            selectingNewReference = true;
        elseif strcmp(fh.SelectionType, 'alt')
            % Do nothing on right click
        elseif strcmp(fh.SelectionType, 'extend') 
            changeNumBins(fh)
        end
        drawSquare(fh)
        selectingNewReference = false;
    end

    function displaySubpatches(im, k)
        imshow(im);
        hold on;
        [X,Y]=meshgrid(0:k:2*r+2);
        plot(X,Y,'k');
        plot(Y,X,'k'); axis off
    end

    function displayScores(numCols, l, dist)
        [X,Y]=meshgrid(0:numCols + 1);
        hold on;
        plot(X,Y,'k');
        plot(Y,X,'k'); axis off;
        for i = numCols:-1:0
            for j = 0:numCols
                text(i, j+0.5, num2str(1 - sum(dist(l:l+C))/(C+1), '%.2f'));
                l = l + 4;
            end
        end
    end

    function enc = subpatchesEncode(x,y)
        l = 1;
        for k = 2.^(0:2)
            for i = (-k+1):2:k
                for j = (-k+1):2:k
                    Dsub = abs(xx-(x - i*ceil(r/k))) <= ceil(r/k) & abs(yy-(y - j*ceil(r/k))) <= ceil(r/k);
                    imgSubpatch = imgLab(Dsub,:);
                    textureSubpatch = imgTexture(Dsub);
                    imgSubpatch = cat(2, imgSubpatch, textureSubpatch);
                    enc((C+1)*l-C:(C+1)*l, :) = patchEncoding(binImage(imgSubpatch,numBins),'hist-normalized',numBins);
                    l = l + 1;
                end
            end
        end
    end

    function showSimilarityHeatmap(encRefPatch)
        SimMap = zeros(H, W);
        for i = r:2*r:W-r
            for j =r:2*r:H-r
            patchEncode = subpatchesEncode(i, j);
            dist = histogramDistance(encRefPatch, patchEncode, 'chi2');
            SimMap(j-r+1:j+r, i-r+1:i+r) = (1 - (w1*dist(1:4) + w2*dist(5:20) + w3*dist(21:84)))*ones(2*r,2*r);
            end
        end
        figure(2)
        heatmap(SimMap)
        title('Similarity to reference patch heat map')
    end


    function showHomogeneityHeatmap()
        HomMap = zeros(H, W);
        for i = r:2*r:W-r
            for j =r:2*r:H-r
            patchEncode = subpatchesEncode(i, j);
            h1 = repmat(patchEncode(1:4),1, 16);
            hom = histogramDistance(h1, patchEncode(21:84), 'chi2');
            HomMap(j-r+1:j+r, i-r+1:i+r) = (1 - max(hom))*ones(2*r,2*r);
            end
        end
        figure(3)
        heatmap(HomMap)
        title('Homogeneity heat map')
    end

end
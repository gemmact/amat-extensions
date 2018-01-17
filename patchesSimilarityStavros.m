function fh = patchesSimilarityStavros(imgRGB)
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
fh = figure(1); imshow(imgRGB);
set(fh, 'WindowButtonMotionFcn', @changePoint);
set(fh, 'WindowButtonDownFcn',   @changeReferencePatch);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
[H,W,C] = size(imgRGB);
tmap    = textonMap(imgRGB, numBins); 
tmap    = reshape(tmap,H*W,[]);
imgRGB  = reshape(imgRGB, [], C);
imgLab  = rgb2labNormalized(imgRGB);
refPatch= zeros(2*r+1);
refxy   = [];
isReferenceFrameSelected = false;
[xx,yy] = meshgrid(1:W,1:H);


    function drawSquare(fh)
        % Get point coordinates and check for validity
        x = round(fh.CurrentAxes.CurrentPoint(1,1));
        y = round(fh.CurrentAxes.CurrentPoint(1,2));
        if x < 1 || x > W || y < 1 || y > H
            title('You are outside of the figure'); drawnow; return
        end
        if x-r < 1 || x+r > W || y-r < 1 || y+r > H
            title('Disk crosses the image boundary'); drawnow; return
        end
        
        % Disk logical mask
        D = abs(xx-x) <= r & abs(yy-y) <= r;

        % The dssim metric should be used on RGB data
        if strcmp(errorType, 'dssim')
            imgPatch = imgRGB(D,:);
        else
            imgPatch = imgLab(D,:);
        end
        
        if isReferenceFrameSelected
            % Keep plotting the reference frame
            rectangle('Position', [xref-r, yref-r, r, r], 'EdgeColor', 'g');
            % Encode the patch and compute error
            encPatch = patchEncoding(binImage(imgPatch,numBins),'hist-normalized',numBins);            
            %err = histogramDistance(refPatch, imgPatch)

            % Disable annoying docking error that clutters the command line
            if strcmp(fh.WindowStyle, 'docked')
                warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
            end
            % Display image and then restore original patch
            % @Gemma: Your visualization code
            disp('I am visualizing')
            drawnow;
        else
            % Simply plot square frame of the patch
            refPatch = imgPatch;
            refxy = [x, y];
            rectangle('Position', [x-r, y-r, r, r], 'EdgeColor', 'g');
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
        elseif strcmp(fh.SelectionType, 'alt')
            % Do nothing on right click
        elseif strcmp(fh.SelectionType, 'extend') 
            changeNumBins(fh)
        end
        drawSquare(fh)
    end
end
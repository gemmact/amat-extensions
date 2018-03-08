function fh = patchesSimilarity(imgRGB, method)
% PATCHESSIMILARITY Choose a reference patch and compare it to any other
%   patch of the same size in the image. It also displays the homogeneity heatmap of the image
%   and the similarity heatmap to the reference patch. This tool allows the 
%   user to interactively change various parameters, such as the size of
%   the patch and whether its shaped as a disk or a square. The supported functions include:
% 
%   [hover mouse over figure]: change the coordinates of the disk/square center.
%   [left click]: select the square/disk centered at the mouse coordinates as the reference patch.
%   [scroll wheel]: up/down (decrease/increase) the radius of the disk/square.
%
%
%   methods supported:
%   'disk'
%   'square'
%
%   See also: imageEncoding


if nargin < 2, method = 'disk'; end

% Default parameters
r = 20;
numBins = 32; % used for histogram encodings
methods.error    = {'se','mse','nmse','rse','rmse','nrmse','dssim'};
methods.encoding = {'average','hist-smirnov','hist-expectation','hist-mode'};
encodingType  = 'hist-normalized';
errorType = 'se';
errorCounter = find(strcmp(methods.error, errorType));
encodingCounter = find(strcmp(methods.encoding, encodingType));

% Plot figure and set callbacks
fh = figure(1);
set(fh, 'Position', [100 100 1100 500])
subplot(3,6,[1 15]);
imshow(imgRGB);
im = imgRGB;
set(fh, 'WindowButtonMotionFcn', @changePoint);
set(fh, 'WindowButtonDownFcn',   @changeReferencePatch);
set(fh, 'WindowScrollWheelFcn',  @changeRadius);
[H,W,C] = size(imgRGB);
imgGray = rgb2gray(imgRGB); 

imgTexture = textonMap(imgGray, numBins);
figure(5); imagesc(imgTexture); axis off image;
%imgTexture    = reshape(imgTexture,H*W,[]);
%tmap    = textonMap(imgRGB, numBins); 
%tmap    = reshape(tmap,H*W,[]);
imgRGB  = reshape(imgRGB, [], C);
imgLab  = rgb2labNormalized(imgRGB);
%imgLabMatrix = reshape(imgLab, [H,W,C]);
refPatch= zeros(2*r+1);
refxy   = [];
encPatch = zeros((C+1)*(1 + r + 1 + (r + 1)/4), numBins);
encRefPatch = zeros((C+1)*(1 + r + 1 + (r + 1)/4), numBins);
isReferenceFrameSelected = false;
selectingNewReference = false;
changingRadius = false;
[xx,yy] = meshgrid(1:W,1:H);
w1 = 0.05*ones(1,4);
w2 = 0.01875*ones(1,16);
w3 = 0.0078125*ones(1,64);


%Show a heatmap with the homogeneity values of all the patches in the image
encImg = imgEncode(method);
%encImgDisk = imgEncodeDisk();
showHomogeneityHeatmap();


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
            switch method
                case 'square'
                    rectangle('Position', [refxy(1)-r, refxy(2)-r, 2*r, 2*r], 'EdgeColor', 'g');
                case 'disk'
                    drawDisk(refxy(1), refxy(2), r, 'g');
            end
            % Encode patch and subpatches
            %%encPatch = subpatchesEncode(y, x, 3);
            %encPatch = subpatchesEncode3(y, x, 1, 1, 3);
            encPatch = subpatchesEncode(y, x, 1, 1, 4, 1:4);

            % Disable annoying docking error that clutters the command line
            if strcmp(fh.WindowStyle, 'docked')
                warning('off','images:imshow:magnificationMustBeFitForDockedFigure')
            end
            
            
            % Display image and compare/save patches
            if (selectingNewReference || changingRadius)
                
                if(changingRadius)
                    encImg = imgEncode();
                    Dref = abs(xx-refxy(1)) <= r & abs(yy-refxy(2)) <= r;
                    if strcmp(errorType, 'dssim')
                        imgRefPatch = imgRGB(Dref,:);
                    else
                        imgRefPatch = imgLab(Dref,:);
                    end
                else
                    imgRefPatch = imgPatch;
                end
                
                
                %save reference patch
                refPatch = reshape(imgRefPatch, [2*r + 1, 2*r + 1, 3]);  
                
                % Encode patch and subpatches
                %%encRefPatch = subpatchesEncode(refxy(2), refxy(1), 3);
                %encRefPatch = subpatchesEncode3(refxy(2),refxy(1),1,1, 3);
                encRefPatch = subpatchesEncode(refxy(2), refxy(1), 1, 1, 4, 1:4);
                
                %Calculate homogeneity score of reference patch
                homRefScore = computeHomogeneity(refxy(1), refxy(2));
                
                %plot patch of reference and subpatches to study
                figure(1)
                
                subplot(3,6,4, 'replace');
                displaySubpatches(refPatch, 0);
                title('Reference patch');
                
                subplot(3,6,10, 'replace');
                displaySubpatches(refPatch, r + 1);
                
                subplot(3,6,16, 'replace');
                displaySubpatches(refPatch, (r+1)/2);
                h = xlabel('');     pos = get(h,'Position'); delete(h)
                h = title(['Hom. score: ' num2str(homRefScore)]); set(h,'Position',pos);
                
                
                
                
            end 
            if (changingRadius || ~selectingNewReference)
                %save patch to compare with reference one
                comparePatch = reshape(imgPatch, [2*r + 1, 2*r + 1, 3]);
 
                %plot square we are comparing to reference one
                switch method
                    case 'square'
                        rectangle('Position', [x-r, y-r, 2*r, 2*r], 'EdgeColor', 'b');
                    case 'disk'
                        drawDisk(x, y, r, 'b');
                end
                
                %Compute homogeneity score of patch
                homScore = computeHomogeneity(x, y);

                
                %plot patch and subpatches to study
                
                subplot(3,6,6, 'replace')
                displaySubpatches(comparePatch, 0);
                title('Second patch');
                    
                subplot(3,6,12, 'replace')
                displaySubpatches(comparePatch, r + 1);

                subplot(3,6,18, 'replace')
                displaySubpatches(comparePatch, (r+1)/2);
                h = xlabel('');     pos = get(h,'Position'); delete(h)
                h = title(['Hom. score: ' num2str(homScore)]); set(h,'Position',pos);
                

                
                %compare patches
                dist = histogramDistance(encRefPatch, encPatch, 'chi2');
                
                %Calculate weighted similarity between the two patches        

                simScore = (1 - (w1*dist(1:4) + w2*dist([5:8 25:28 45:48 65:68]) + w3*dist([9:24 29:44 49:64 69:84])));
                
                %plot similarity scores between histograms
                
                subplot(3,6,5, 'replace')
                displayScores(0, 1, dist, 0);
                

                subplot(3,6,11, 'replace')
                displayScores(1, 5, dist, 20);

                subplot(3,6,17, 'replace')
                displayScores(3, 9, dist, 4);
                h = xlabel('');     pos = get(h,'Position'); delete(h)
                h = title(['Similarity score: ' num2str(simScore)]); set(h,'Position',pos);
                

            end
            
            
            drawnow;
        else
            % Simply plot square frame of the patch
            switch method
                case 'square'
                    rectangle('Position', [x-r, y-r, 2*r, 2*r], 'EdgeColor', 'b');
                case 'disk'
                    drawDisk(x, y, r, 'b');
            end
        end        
    end


    function changePoint(fh,~)
        drawSquare(fh);
    end

    function changeRadius(fh,callbackData)
        changingRadius = true;
        r = min(min(H,W)/2, max(1, r + callbackData.VerticalScrollCount));
        showHomogeneityHeatmap()
        drawSquare(fh);
        changingRadius = false;
        %Create heatmap of similarity values between the reference
        %patch and the other ones      
        showSimilarityHeatmap(encRefPatch);
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
        
        %Create heatmap of similarity values between the reference
        %patch and the other ones     
        showSimilarityHeatmap(encRefPatch);
    end

    function displaySubpatches(im, k)
        switch method
            case 'disk'
                imshow(im.*disk(r));
            case 'square'
                imshow(im);
        end
        hold on;
        if (k>0)
            [X,Y]=meshgrid(0:k:2*r+2);
            plot(X,Y,'k');
            plot(Y,X,'k'); axis off
        end
    end

    function displayScores(numCols, l, dist, k)
        [X,Y]=meshgrid(0:numCols + 1);
        hold on;
        plot(X,Y,'k');
        plot(Y,X,'k'); axis off;
        for j = numCols:-1:0        
            for i = 0:numCols
                text(i, j+0.5, num2str(1 - sum(dist(l:l+C))/(C+1), '%.2f'));
                %l = l + 4;
                l = l+ k;
                if(i == (numCols) && mod(j,2) == 1 && numCols>=3) 
                    l = l - (2*k);
                end
                if (~(i == (numCols-1) && mod(j,2) == 0))
                    if (k == 4) 
                        k = 16;
                    else if (k == 16) 
                            k = 4;
                        end
                         if (k == 8)
                             k = 4;
                        end
                    end
                else
                    k = 8;
                end
            end
        end
    end


    function encod = subpatchesEncode(x,y,k,v, k_end, chan)
       encod1 = reshape(encImg(x, y, chan, :, v), [length(chan), numBins]);
       if(k<k_end)
           switch method
               case 'disk'
                   encod2 = subpatchesEncode(x - ceil(r/(2*k)),y - ceil(r/(2*k)), 2*k, 4*v-2, k_end, chan);
                   encod3 = subpatchesEncode(x - ceil(r/(2*k)),y + ceil(r/(2*k)), 2*k, 4*v-1, k_end, chan);
                   encod4 = subpatchesEncode(x + ceil(r/(2*k)),y - ceil(r/(2*k)), 2*k, 4*v, k_end, chan);
                   encod5 = subpatchesEncode(x + ceil(r/(2*k)),y + ceil(r/(2*k)), 2*k, 4*v+1, k_end, chan);
               case 'square'
                   encod2 = subpatchesEncode(x - ceil(r/(2*k)),y - ceil(r/(2*k)), 2*k, v+1, k_end, chan);
                   encod3 = subpatchesEncode(x - ceil(r/(2*k)),y + ceil(r/(2*k)), 2*k, v+1, k_end, chan);
                   encod4 = subpatchesEncode(x + ceil(r/(2*k)),y - ceil(r/(2*k)), 2*k, v+1, k_end, chan);
                   encod5 = subpatchesEncode(x + ceil(r/(2*k)),y + ceil(r/(2*k)), 2*k, v+1, k_end, chan);  
           end
           encod = [encod1; encod2; encod3; encod4; encod5];
       else
           encod = encod1;
       end
    end

    function encod = subpatchesEncodeDiag(x,y,k,v, v_end, chan, i)
       encod1 = reshape(encImg(x, y, chan, :, v), [length(chan), numBins]);
       if(v<v_end)
           switch method
               case 'disk'
                   encod2 = subpatchesEncodeDiag(x + i*ceil(r/(4*k)),y - ceil(r/(2*k)), 2*k, v+1, v_end, chan, i);
                   encod3 = subpatchesEncodeDiag(x - ceil(r/(2*k)),y - i*ceil(r/(4*k)), 2*k, v+2, v_end, chan, i);
                   encod4 = subpatchesEncodeDiag(x + ceil(r/(2*k)),y + i*ceil(r/(4*k)), 2*k, v+3, v_end, chan, i);
                   encod5 = subpatchesEncodeDiag(x - i*ceil(r/(4*k)),y + ceil(r/(2*k)), 2*k, v+4, v_end, chan, i); 
           end
           encod = [encod1; encod2; encod3; encod4; encod5];
       else
           encod = encod1;
       end
    end


    function encod = subpatchesEncode3(x,y,k,v, v_end)
       encod1 = reshape(encImg(x, y, 1:4, :, v), [C+1, numBins]);
       if(v<v_end)
           encod2 = subpatchesEncode3(x - ceil(r/(2*k)),y - ceil(r/(2*k)), 2*k, v+1, v_end);
           encod3 = subpatchesEncode3(x - ceil(r/(2*k)),y + ceil(r/(2*k)), 2*k, v+1, v_end);
           encod4 = subpatchesEncode3(x + ceil(r/(2*k)),y - ceil(r/(2*k)), 2*k, v+1, v_end);
           encod5 = subpatchesEncode3(x + ceil(r/(2*k)),y + ceil(r/(2*k)), 2*k, v+1, v_end);
           encod = [encod1; encod2; encod3; encod4; encod5];
       else
           encod = encod1;
       end
    end



    function showSimilarityHeatmap(encRefPatch)
        SimMap = zeros(H, W);
        
        %Encode all patches of size 2r x 2r and their subpatches and
        %calculate similarity with reference patch
        for i = r+1:W-r-1
            for j =r+1:H-r-1    
            %%patchEncode = subpatchesEncode(j, i, 3);
            patchEncode = subpatchesEncode(j, i, 1, 1, 4, 1:4);
            dist = histogramDistance(encRefPatch, patchEncode, 'chi2');
            %SimMap(j, i) = (1 - (w1*dist(1:4) + w2*dist(5:20) + w3*dist(21:84)));
            SimMap(j, i) = (1 - (w1*dist(1:4) + w2*dist([5:8 25:28 45:48 65:68]) + w3*dist([9:24 29:44 49:64 69:84])));
            end
        end
        
        %Display heatmap with similarity values to reference patch
        figure(3)
        heatmap(SimMap, 'GridVisible', 'off', 'XDisplayLabels', strings(W,1), 'YDisplayLabels', strings(H,1))
        title('Similarity to reference patch heat map')
    end

    function homScore = computeHomogeneity(i, j)
            patchEncode = subpatchesEncode(j, i, 1, 1, 2, 1:4);
            %patchEncodeTex = subpatchesEncode(j, i, 1, 1, 2,4);
            h1 = repmat(patchEncode(1:4,:),4, 1);
            %h2 = repmat(patchEncodeTex(1,:),4, 1);
            %hom = histogramDistance(h1, patchEncode(4:15, :), 'chi2');
            %hom2 = histogramDistance(h2, patchEncodeTex(2:5, :), 'chi2');
            %HomMap(j, i) = (1 - max([hom;hom2]));
            %patchEncode = subpatchesEncode3(j, i, 1, 1, 2);
            %h1 = repmat(patchEncode(1:4,:),4, 1);
            hom = histogramDistance(h1, patchEncode(5:20, :), 'chi2');
            hom = sum(0.25*reshape(hom, [4, 4]));
            switch method
                case 'disk'
                for k = -1:1   
                    patchEncodeDiag = subpatchesEncodeDiag(j, i, 1, 22 + (k+1)*5, 23 + (k+1)*5, 1:4, k);
                    h1Diag = repmat(patchEncodeDiag(1:4,:),4, 1);
                    homDiag = histogramDistance(h1Diag, patchEncodeDiag(5:20, :), 'chi2');
                    homDiag = sum(0.25*reshape(homDiag, [4, 4]));
                    %homTot = [(hom+homDiag)/2 , (hom+[homDiag(2) homDiag(4) homDiag(1) homDiag(3)])/2];
                    homAux = [homDiag(2) homDiag(4) homDiag(1) homDiag(3)];
                    homTot = [max(hom(1:4), homDiag), max(hom(1:4), homAux)];
                    l = length(hom);
                    if (l > 4)
                        for t = 2:(l/4)
                        homTot = [homTot, max(hom(4*(t-1)+1:4*t), homAux)];
                        end
                    end
                    hom = homTot;
                end
                
                case 'square'
                homTot = hom;
            end
            
            homScore = (1 - max(homTot));
    end


    function showHomogeneityHeatmap()
        HomMap = zeros(H, W);
        
        %Encode histograms for all 2rx2r patches in the image and its 
        %subpatches and compare the finest with the coarsest scales to
        %compute homogeneity scores
        for i = r+1:W-r-1
            for j =r+1:H-r-1
            HomMap(j, i) = computeHomogeneity(i, j);
            end
        end
        
        %Display homogeneity heat map
        figure(2)
        heatmap(HomMap, 'GridVisible', 'off', 'XDisplayLabels', strings(W,1), 'YDisplayLabels', strings(H,1))
        title('Homogeneity heat map')
    end

    function enc = imgEncode(method)
        switch method
            case 'square'
               Dsub = cell(3,1);
               %DsubTex = cell(3,1);
               for i = 1:3 
                kk = 2^(i-1);
                Dsub{i} = ones(2*ceil(r/kk));
                %DsubTex{i} = ones(2*ceil(r/kk));
               end
            case 'disk'
                Dsub = cell(36,1);
                %Dsubdiag = cell(5, 1);
                %DsubTex = cell(5,1);
                Dsub{1}=disk(r);
                Dsub{2}=disk(r, 'quadrant4');
                Dsub{3}=disk(r, 'quadrant3');
                Dsub{4}=disk(r, 'quadrant2');
                Dsub{5}=disk(r, 'quadrant1');
                for i = 2:5
                    Dsub{1+ 4*(i-1) + 1} = Dsub{i}(floor(r/2):r, floor(r/2):r);
                    Dsub{1+ 4*(i-1) + 2} = Dsub{i}(floor(r/2):r, 1:ceil(r/2));
                    Dsub{1+ 4*(i-1) + 3} = Dsub{i}(1:ceil(r/2), floor(r/2):r);
                    Dsub{1+ 4*(i-1) + 4} = Dsub{i}(1:ceil(r/2), 1:ceil(r/2));
                    %Dsub{21+i}=imrotate(T, -90*i);
                end
                
                for i = 1:3
                    %U = triu(imrotate(Dsub{1}, (i-2)*22.5));
                    U = triu(Dsub{1});
                    L = imrotate(U, 90);
                    T = L.*U;
                    T = imrotate(T, (i-2)*22.5);
                    [A,B] = size(T);
                    T = T(1:ceil(A/2), :);
                    
                    Dsub{22 + (i-1)*5}=disk(r);
                    Dsub{23 + (i-1)*5} = imrotate(T, 270);
                    Dsub{24 + (i-1)*5} = imrotate(T, 180);
                    Dsub{25 + (i-1)*5} = T;
                    Dsub{26 + (i-1)*5} = imrotate(T, 90);
                end
                

            otherwise
                error('Method not supported')
        end
       %imgSubpatch = binImage(reshape(imgLab, [H,W,C]), numBins);
       %enc = cat(3, imageEncoding(imgSubpatch, Dsub, 'hist-normalized', numBins),imageEncoding(imgTexture, DsubTex, 'hist-normalized', numBins));
        
       imgSubpatch = cat(3, binImage(reshape(imgLab, [H,W,C]), numBins), imgTexture);
       size(Dsub)
       enc = imageEncoding(imgSubpatch, Dsub, 'hist-normalized', numBins);
       %enc = cat(3, imageEncoding(imgSubpatch, Dsub, 'hist-normalized', numBins),imageEncoding(imgSubpatch, Dsubdiag, 'hist-normalized', numBins));
    end



    function drawDisk(xCenter, yCenter, r, col)
        theta = 0 : 0.01 : 2*pi;
        x = r * cos(theta) + xCenter;
        y = r * sin(theta) + yCenter;
        hold on
        plot(x, y, col);
        hold off
    end

end
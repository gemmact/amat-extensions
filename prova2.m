%% Test black and white patterns
B = 32;
sz = [100,100];
black = zeros(sz);
white = ones(sz);
half = [black(:,1:sz(2)/2), white(:,1:sz(2)/2)];
checker = white; checker(1:2:end-1,1:2:end) = 0; checker(2:2:end,2:2:end) = 0;
checkerinv = 1-checker;
gray = mean(half(:))*ones(sz);
hblack = patchEncoding(binImage(black(:),B),'hist-normalized',B);
hwhite = patchEncoding(binImage(white(:),B),'hist-normalized',B);
hhalf = patchEncoding(binImage(half(:),B),'hist-normalized',B);
hchecker = patchEncoding(binImage(checker(:),B),'hist-normalized',B);
hcheckerinv = patchEncoding(binImage(checkerinv(:),B),'hist-normalized',B);

thalf = textonMap(half,B);
tchecker = textonMap(checker,B);
htexhalf = patchEncoding(thalf(:),'hist-normalized',B);
htexchecker = patchEncoding(tchecker(:),'hist-normalized',B);

%% Test different parameters for extracting textons
B = 32;
tmap1 = textonMap(half,B,'global',6,1,1,sqrt(2),2);
tmap2 = textonMap(checker,B,'global',6,1,2,sqrt(2),2);
% tmap1 = computeTextons(fbRun(fbCreate(6,1,1,2),half),B);
% tmap2 = computeTextons(fbRun(fbCreate(6,1,1,2),checker),B);
h1 = patchEncoding(tmap1(:),'hist-normalized',B);
h2 = patchEncoding(tmap2(:),'hist-normalized',B);
figure(1); bar(h1); figure(2); bar(h2)
figure(3); imagesc(tmap1); axis off image;
figure(4); imagesc(tmap2); axis off image;
% figure(3); imshow(imgRGB);

% textons = unitex(fbCreate(6,1,1,3),32);

%% Test different pre-processing for input image.
imggaussf = imgaussfilt(imgRGB,1);
imgmedian = zeros(size(imgRGB)); 
for c=1:size(imgRGB,3), imgmedian(:,:,c) = medfilt2(imgRGB(:,:,c),[5 5]); end
imgavg = imfilter(imgRGB,fspecial('disk',3));
figure;
subplot(221); imshow(imgRGB); title('Original')
subplot(222); imshow(imggaussf); title('Gaussian')
subplot(223); imshow(imgmedian); title('Median')
subplot(224); imshow(imgavg); title('Average (disk)')



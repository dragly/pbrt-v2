% read the image
content = dlmread('sprite1.exr.vsd');

% get the image resolution 
imageLength = size(content);
xResolution = content(imageLength(1), 1);
yResolution = content(imageLength(1), 2);
pro
% get the image 
vsdImage = zeros(xResolution, yResolution);

% fill the data in the image 
for i = 1:imageLength(1)
   x = content(i, 1) + 1;
   y = content(i, 2) + 1;
   value = content(i, 3);
   vsdImage(x, y) = value;   
end

figure;
imagesc(vsdImage);

%function [im1,im2,im3,im4,centers, centroid,indexOfRemainingRegions] = fcmDemo3()
V = niftiread('patient001_frame01_gt.nii.gz');
c = V(:,:,9);
J = imrotate(c,-90,'bilinear','crop');
I = flip(J ,2);
imshow(I)
n_clusters = 2;
imData= reshape(I,[],1);
imData= double(imData);

[centers,U,obj]= fcm(imData,n_clusters);

c1 = U(1,:);

c2= U(2,:);

im1 = reshape(c1,size(I));
im2 = reshape(c2,size(I));
im3 = im2bw(I, graythresh(I));
im4 = imfill(im3, 'holes');

%subplot(2,2,1),imshow(im1),title('Cluster 1');
%subplot(2,2,2),imshow(im2),title('Cluster 2'); 
%subplot(2,2,3),imshow(im3),title('Cluster 3');
%subplot(2,2,4),imshow(im4),title('Cluster 4'); 

originalImage = im2; 
binaryImage= im4;
captionFontSize = 14;

% Label each blob so we can make measurements of it
labeledImage = bwlabel(im4, 8); 

% Let's assign each blob a different color to visually show the user the distinct blobs.
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
% coloredLabels is an RGB image.  We could have applied a colormap instead (but only with R2014b and later)
%subplot(4, 4, 1);
 figure;
imshow(I,[]);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
%caption = sprintf('(a) Original Image');
title(caption, 'FontSize', captionFontSize);
% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.

% (b) Partitioned Regions
%subplot(4, 4, 2);
 figure;
imshow(binaryImage,[]);
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
caption = sprintf('(b) Partitioned Regions');
title(caption, 'FontSize', captionFontSize);
% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.


blobMeasurements = regionprops(labeledImage, originalImage, 'all');
numberOfBlobs = size(blobMeasurements, 1);
% bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.

 figure;
imshow(I,[]);
title('(c) Delineated Contours', 'FontSize', captionFontSize); 
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
end
hold off;

textFontSize = 2;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
blobECD = zeros(1, numberOfBlobs);
% Print header line in the command window.
fprintf(1,'Blob #      Mean Intensity  Area   Perimeter    Centroid       Diameter\n');
% Loop over all blobs printing their measurements to the command window.
area=zeros(numberOfBlobs,1);%Global variable to use for later 
centroid=zeros(numberOfBlobs,2);
circularity=zeros(numberOfBlobs,1);
blobBox=zeros(numberOfBlobs,4);

for k = 1 : numberOfBlobs           % Loop through all blobs.
	% Find the mean of each blob.  (R2008a has a better way where you can pass the original image
	% directly into regionprops.  The way below works for all versions including earlier versions.)
	thisBlobsPixels = blobMeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
	meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
	meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a
	
	blobArea = blobMeasurements(k).Area;		% Get area.
	blobPerimeter = blobMeasurements(k).Perimeter;		% Get perimeter.
	blobCentroid = blobMeasurements(k).Centroid;		% Get centroid one at a time

	blobECD(k) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
	fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD(k));
	% Put the "blob number" labels on the "boundaries" grayscale image.
	%text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold');
    % For Step 2.1c
    area(k) = blobArea;
    centroid(k,:) = blobCentroid;
    circularity(k)=(4*pi*(blobArea/blobPerimeter^2));
    blobBox(k,:) = blobMeasurements(k).BoundingBox;
end

% Step 2.1d
Area_Tol= 374;
Round_Tol= 0.6;
Dist_Tol=  27 ;
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.

 figure;
imshow(I,[]);
title('(d) After Area', 'FontSize', captionFontSize); 
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
    if (area(k)>=Area_Tol)
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
    end
end
hold off;

% Step 2.1e 
% Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.

 figure;
imshow(I,[]);
title('(e) After Circularity Filter', 'FontSize', captionFontSize); 
axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);
numberOfRemainingRegions=0;
for k = 1 : numberOfBoundaries
    if (area(k)>=Area_Tol)
        if(circularity(k) >=Round_Tol)
            thisBoundary = boundaries{k};
            plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
            plot(centroid(k,1),centroid(k,2),'b*')
            numberOfRemainingRegions = numberOfRemainingRegions+1;
        end
    end
end
hold off;

% Step 2.1f 
indexOfRemainingRegions=zeros(numberOfRemainingRegions);
i=1;
for k = 1 : numberOfBoundaries % get the index Of Remaining Regions
 if (area(k)>=Area_Tol)
        if(circularity(k) >=Round_Tol)
            indexOfRemainingRegions(i)=k;
            i=i+1;
        end
 end
end
xyCentroid =0;
for k = 1 : i-1
    xyCentroid(k,1) = centroid(indexOfRemainingRegions(k,1),1);
    xyCentroid(k,2) = centroid(indexOfRemainingRegions(k,1),2);
end

% Calculate the distance from each region centroid to the image borders
if(xyCentroid ~=0)
    xyImageEdge = size(I);
     distancesR = pdist2(xyCentroid, xyImageEdge);
     distancesL = pdist2(xyCentroid, [xyImageEdge(1,1),1]);
     distances = distancesR+distancesL;
     % Choose the farthest region from the image border 
     HighestDistance = max(distances);
     HighestDistanceIndex=0;
    for k = 1 : i-1
      if(HighestDistance == distances(k))
          HighestDistanceIndex=k;
      end
    end
    % Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.

     figure;
    imshow(I,[]);
    %title('(d) Frame 2', 'FontSize', captionFontSize); 
    axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
    hold on;
    boundaries = bwboundaries(binaryImage);
    numberOfBoundaries = size(boundaries, 1);
    numberOfRemainingRegions=0;
    for k = 1 : numberOfBoundaries
        if (indexOfRemainingRegions(HighestDistanceIndex)==k)
                thisBoundary = boundaries{k};
                plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
                plot(centroid(k,1),centroid(k,2),'b*')
                numberOfRemainingRegions = numberOfRemainingRegions+1;
        end
    end
    hold off;
 
end 

% Step 2.1g

 figure;
imshow(I,[]);
    %title('(b) ROI', 'FontSize', captionFontSize); 
    axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
    hold on;
    boundaries = bwboundaries(binaryImage);
    numberOfBoundaries = size(boundaries, 1);
    for k = 1 : numberOfBoundaries
        if (indexOfRemainingRegions(HighestDistanceIndex)==k)
                plot(centroid(k,1),centroid(k,2),'g*')
                rectangle('Position',[blobBox(k,1),blobBox(k,2),blobBox(k,3),blobBox(k,4)],'EdgeColor','g','LineWidth',1)
                rectangle('Position',[blobBox(k,1)/1.14,blobBox(k,2)/1.23,blobBox(k,3)*1.6,blobBox(k,4)*1.6],'EdgeColor','r','LineWidth',1)
        end
    end
    hold off;

%end

     %Function to calculate approximate size of clusters using an intensity
%cutoff. This code necessarily will have to get more sophisticated as we
%start to more carefully analyze these images. This code is not intended to
%analyze data at a single cell level, but rather give bulk features of the
%aggregates.
%
%% 
% <intensitySegmentLargeClusters: Analysis of aggregate formation for a time series>
%     Copyright (C) 2018, Matthew Jemielita
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%     

function [lm,cc,cm,frac] = intensitySegmentLargeClusters(im, th, minObjSize, minObjVol,mag)
%Contains different subroutines for analyzing formation of clusters in LCD
%or HCD aggregates.

lm = [];cc = []; cm= [];
frac = [];
switch mag
    case '10x'
        [lm,cc,cm] = segment_10x(im, th, minObjSize, minObjVol); 
        
    case '63x'
        %Intensity-based segmentation of LCD-locked clusters
        lm = segment_63x(im, th, minObjSize, minObjVol);
        
    case '10x volume'
        [lm,cc] = segmentVolume_10x(im, th, minObjSize, minObjVol);

end

end

function [lm,cc,cm] = segment_10x(im, th, minObjSize, minObjVol)
%Manually remove regions that are outside the well. Also useful for
%preventing false positives from the ring of intenstiy around the outside
%of the well.
figure; imshow(max(im,[],3),[]);
h = imellipse;
pos = wait(h);

mask = createMask(h);

%Coarse segmentation
bw = im>th;

%Apply mask to remove things touching the boundary;
bw = bw.*repmat(mask, [1,1, size(bw,3)]);
bw = bw>0;%Make logical again, otherwise we get weird bugs.
%Cleaning up image to remove noise, etc.
fprintf(1, 'Cleaning up image');

%% Fill in holes and clean up objects sticking to the border
for i= 1:size(im,3)
    bw(:,:,i) = bwareaopen(bw(:,:,i),minObjSize);
    bw(:,:,i) = imclearborder(or(bw(:,:,i),~mask));
    bw(:,:,i) = imfill(bw(:,:,i),'holes');
    
  %  bw(:,:,i) = imclose(bw(:,:,i), strel('disk', 1));
    fprintf(1, '.');
end
fprintf(1, 'done!\n');

%%
cc = bwconncomp(bw, 18);
cc = regionprops(cc, 'Area', 'Centroid', 'PixelIdxList' ,'BoundingBox');
%%
zSize = [cc.BoundingBox];
zSize = reshape(zSize,6,length(zSize)/6);
zSize = zSize(6,:);

ind = zSize>1;
cc = cc(ind);

ind = [cc.Area]>minObjVol;

cc = cc(ind);

%% Make a 3-d label matrix for the segmented images
lm = zeros(size(bw));

for i=1:length(cc)
lm(cc(i).PixelIdxList) = i;

end
%%
cm = lines(length(cc));

end

function [finalLinkage, bwlabelout] = segmentVolume_10x(im, th, minObjSize, minObjVol)
%Segmentation of aggregates seen in a HCD locked background.
%minObjVol coded in but not currently used
%Coarse segmentation
bw = im>th;

%% Cleaning up image to remove noise, etc.
fprintf(1, 'Cleaning up image');
bw = medfilt3(bw);
%% Fill in holes and clean up objects sticking to the border
for i= 1:size(im,3)  
    bw(:,:,i) = bwareaopen(bw(:,:,i),minObjSize,4);
    bw(:,:,i) = imfill(bw(:,:,i),'holes');
    bw(:,:,i) = imclearborder(bw(:,:,i));
    bw(:,:,i) = imopen(bw(:,:,i), strel('disk', 4));
    fprintf(1, '.');
end
fprintf(1, 'done!\n');


%% Observation: All aggregates are basically convex in shape. Let's use this
%observation to assess whether we have undersegmented touching aggregates.
bwlabel = zeros(size(bw));

for i=1:size(im,3)
    bwlabelout(:,:,i) = UseConvexProperty(bw(:,:,i));
end


%% Link together 3d structures
 finalLinkage = link3Dvolume(bwlabelout);

end

function bw = segment_63x(im, th, minObjSize, minObjVol)
%Coarse segmentation based on intensity based segmentation, image closure
%and removing small objects
%minObjVol coded in but not currently used
%% Cleaning up image to remove noise, etc.
fprintf(1, 'Cleaning up image');
%% Fill in holes and clean up objects sticking to the border
for i= 1:size(im,3)  
    t = imclose(im(:,:,i), strel('disk', 5));
    t = t > th;
    
    t = ~bwareaopen(~t,1000);
    t = bwareaopen(t, minObjSize,4);
    t = imclearborder(t);
    
    bw(:,:,i)= t;
    fprintf(1, '.');
end
fprintf(1, 'done!\n');
bw = bw>0;
end

function [mask1, mask2] = getMinimalCut(cc_test,mask,i)
%Compute the cut across undersegmented objects that is the minimal
%distance. This is a reliable way to split apart aggregates that are
%touching.
mask(cc_test(i).PixelIdxList) = 1;

b = bwboundaries(mask); b = b{1};
x = b(:,1);y = b(:,2);

%Let's sample approximately every hundred points
indrs = linspace(1,length(x),100);
indrs = round(indrs);
indrs = unique(indrs);
x = x(indrs);
y = y(indrs);

validLine = cell(1,1);

%Assemble all bisections through the region that have a line of sight
%exclusively through the polygon.
% figure; plot(x,y);hold on;
n=1;
for i=1:length(indrs)
    for j=i+1:length(indrs)
        xl = linspace(x(i), x(j));
        yl = y(i) + (y(j)-y(i))/(x(j)-x(i))  *(xl-x(i));
        
        in = inpolygon(xl,yl, x,y);
        in = prod(in);
        
        if(in==1)
            validLine{n}.line = [xl;yl];
            validLine{n}.i = i;
            validLine{n}.j = j;
            validLine{n}.length= sqrt((yl(end)-yl(1))^2  +  (xl(end)-xl(1))^2);
            n = n+1;
        end
        
    end
    fprintf(1, '.');
end

%% Find all of the lines that are realy short
for i=1:length(validLine)
    t(i) = validLine{i}.length;
    d(i) = validLine{i}.i-validLine{i}.j;
    d(i) = sqrt(d(i)^2);
    
end
%% Define the cut as the one that is the minimum length for the middle 50% of index-index distances

ind = (d>25).*(d<75);
ind = logical(ind);

[~,indMin] = min(t(ind));

indVal = find(ind==1);
%Note: This gives a reasonable value for our segmented piece.
indVal = indVal(indMin);
%% For each of these lines, construct the two new polygons that would result if
%we cut the object at those points

ind_i = validLine{indVal}.i;
ind_j = validLine{indVal}.j;

poly1.x = [x(ind_i:ind_j-1); validLine{indVal}.line(1,:)'];
poly1.y = [y(ind_i:ind_j-1); validLine{indVal}.line(2,:)'];

poly2.x = [validLine{indVal}.line(1,:)'; x(ind_j:end); x(1:ind_i)];
poly2.y = [validLine{indVal}.line(2,:)'; y(ind_j:end); y(1:ind_i)];

allpolygon.poly1 = poly1;
allpolygon.poly2 = poly2;


mask1 = poly2mask(allpolygon.poly1.y, allpolygon.poly1.x, size(mask,2), size(mask,1));
mask2 = poly2mask(allpolygon.poly2.y, allpolygon.poly2.x, size(mask,2), size(mask,1));
end

function maskOutput = UseConvexProperty(mask)
        convexThresh = 4000;
        
        cc  = bwconncomp(mask,  4);
        cc = regionprops(cc, 'ConvexArea', 'Area', 'PixelIdxList');
        
        %%
        temp = mask;
        temp(:) = 0; temp = double(temp);
        for i=1:length(cc)
            temp([cc(i).PixelIdxList]) = cc(i).ConvexArea./cc(i).Area;
        end
   
        %% Compute rms difference between convex area and area
        ca_ratio =  sqrt(([cc.Area]-[cc.ConvexArea]).^2 );
        bwRatio =(-1)* zeros(size(temp));
        
        for i=1:length(cc)
            bwRatio([cc(i).PixelIdxList]) = ca_ratio(i);
        end
 
        %Go through objects that have a large convex/non-convex ratio
        %Shrink them, and measure the new ratio. If there's a point where it breaks
        %into two objects
        
        cc_test = cc(ca_ratio> convexThresh);
        
        mask = zeros(size(temp));
        maskFinal = mask;
        n =1;
        for i=1:length(cc_test)     
            [mask1, mask2] = getMinimalCut(cc_test,mask,i);
            maskFinal = maskFinal+ n*mask1;
            maskFinal = maskFinal+ (n+1)*mask2;
            n = n+2;
        end
        fprintf(1, '\n');
        
        
        %% Combining together original clusters and the ones that were cut apart that had high convex ratio
        bwRatio(bwRatio>convexThresh) = 0;
        temp = temp.*(bwRatio>0);
        
        temp2 = temp+maskFinal;
        
        label = unique(temp2(:));
        maskOutput = bwlabel(bwRatio); 
        maxInd = max(maskOutput(:));
        maskFinal = maskFinal + maxInd;
        
        maskOutput = maskOutput + maskFinal;
        maskOutput = imclearborder(maskOutput);
end
    
function finalLinkage = link3Dvolume(bwlabelout)
%%
finalLinkage = zeros(size(bwlabelout));
finalLinkage(:,:,1) = bwlabelout(:,:,1);
linkageInd = 1+ max(finalLinkage(:));
%% Find overlapping objects

for i=2:size(bwlabelout,3)
    
    im = bwlabelout(:,:,i);
    imprev = finalLinkage(:,:,i-1);
    
    %Find centroid of each of these objects-this will be used for matching
    %them from frame to frame
    ccim = regionprops(im, 'Centroid');
    ccprev = regionprops(imprev, 'Centroid');
    
    
    %Find all unique objects in this frame
    ind = unique(im(:));
    ind(ind==0) = [];
    %If frame is empty, skip
    if(isempty(ind))
        continue;
    end
    
    %Construct blank image for this frame that we'll fill with new cluster
    %ID's
    updatedFrame = zeros(size(im));
    
    %Find list of indices that are overlapping with clusters in the previous
    %frame
    overlapBW = (im.*imprev)>0;
    indOverlap = im(overlapBW);
    
    indOverlapPrev = imprev(overlapBW);
    indOverlapUniq = unique(indOverlap);
    
    % For now, link together the first found entry in this list of frames
    % that link together.
    %NOTE: ***************** This is a point where future reworking will
    %probably be necessary to make sure the linking works as well as
    %possible
    %Really what I should do is do this globally using something like the
    %hungarian algorithm used for the assignment problem. Code for this
    %exists on the Danuser lab website-might take a little bit of time to
    %implement into our code.
    
    %Identify index number for clusters that overlap with previous clusters.
    %Update new image with those cluster ID's
    
    possibleLinkedCluster = unique([indOverlap, indOverlapPrev], 'rows');
    
    mappingIndex = zeros(length(indOverlapUniq),2);
    for j=1:length(indOverlapUniq)
        
        thisInd = (possibleLinkedCluster(:,1)==indOverlapUniq(j));
        thisInd = possibleLinkedCluster(thisInd,:);
        
        centrIm = [ccim(thisInd(:,1)).Centroid];
        centrPrev = [ccprev(thisInd(:,2)).Centroid];
        centrIm = reshape(centrIm, 2,length(centrIm)/2);
        centrPrev = reshape(centrPrev, 2, length(centrPrev)/2);
        
        dist = sum(sqrt((centrIm-centrPrev).^2),1);
        [~, closestInd] = min(dist);
        mappingIndex(j,:) = thisInd(closestInd,:);
    end
    
    for j=1:size(mappingIndex,1)
        thisObj = im==mappingIndex(j,1)>0;
        updatedFrame(thisObj) = mappingIndex(j,2);
    end
    
    %Find list of indices that are in clusters that do not overlap with
    %clusters in the previous frame
    
    %linkageInd = linkageInd + max(updatedFrame(:));
    
    indNew = setdiff(ind, indOverlapUniq);
    indUpdated = linkageInd:linkageInd+length(indNew)-1;
    
    for j=1:length(indNew)
        thisObj = im==indNew(j)>0;
        updatedFrame(thisObj) = indUpdated(j);
    end
    
    %Update linkageInd
    linkageInd = linkageInd + length(indNew);
    
    %Update final linked image
    finalLinkage(:,:,i) = updatedFrame;
    fprintf(1, '.');
end
fprintf(1, '\n');
end

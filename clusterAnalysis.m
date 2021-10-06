% <clusterAnalysisTimeSeries: Analysis of aggregate formation for a given data set>
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

classdef clusterAnalysis
   properties
       saveLoc = '';
       saveName = '';
       labelSaveName = '';
       
       genotype = '';
       
       median = '';
       skewness = '';
       
       binsize = 1;
       binlimits = [];
       histdata = [];
       histedges = [];
       %Variables for accessing the raw data that will be segmented
       imageDataLoc = '';
       imageDataInd = 0;
       
       imageDataInd_SecondChannel = 0; %Location of the second image channel
       %Default variables for segmenting 3d objects
       th = 100;
       minObjSize = 100;
       minObjVol = 10000;
       zExt = 2;
       
       labelm = []; %Label matrix of clusters
       
       time = 0;
       
       %Data structure holding info about each of the clusters
       cc
       
   end
   
   methods
       
       function obj = clusterAnalysis(saveLoc, genotype, time)
           %Create instance of class
           obj.saveLoc = saveLoc;
           obj.genotype = genotype;
           
           obj.time = time;
           obj.saveName = [obj.genotype '_t' num2str(obj.time)];
           obj.labelSaveName = [obj.saveName '_labelmatrix.mat'];
       end
       
       function obj = runPipeline(obj)
           %Run pipeline for segmented aggregates
           obj = obj.segmentClusters;
           obj = obj.getRegionProps;
           obj.labelm = [];
           obj.save;
       end
       
       function obj = runPipeline_2color(obj, color_1, color_2)
          %Run pipeline for two color channels, using one as the marker for the extent of the other
          %Use color_1 as the master color to segment the objects.
          %Measure the fractional occupancy in the aggregates found in
          %color_1 channel for the color_2 channel.
          %Pass in the path and filename to either color channel.
          
          %Load color_1 channel image. These are stored as tiff's since
          %bfopen seems to have problems with multi-channel images.
          
          %Preallocate to known size
          im_1 = zeros(2048,2048,51);
          for i=1:51
              im_1(:,:,i) = imread(color_1, 'Index', i);
          end
          im_2 = zeros(2048,2048,51);
          
          for i=1:51
              im_2(:,:,i) = imread(color_2, 'Index', i);
          end
          
          sprintf('Images loaded');
          [obj.labelm, ~, ~,~] = intensitySegmentLargeClusters(im_1,...
              obj.th,obj.minObjSize, obj.minObjVol, '10x volume');
          obj = obj.getRegionProps_2color(im_1, im_2);
          obj.labelm = [];
          obj.save;
          
           
       end
       
       function obj = runPipelineLCD(obj)
           %Run pipeline for segmenting clusters (objects that appear in a
           %LCD-locked background)
           obj = obj.segmentLCDClusters;
           obj = obj.getRegionProps;
           obj.labelm = [];
           obj.save;
       end
       
       function obj = runPipeLineXdsDns(obj)
          %Do we need this?
           obj = obj.segmentClusters_xds_dns;
           obj = obj.getRegionProps;
          obj.labelm = [];
          %obj.constructHistogram;
          obj.save; 
       end
       
       function obj = segmentClusters(obj)
           %Run through the code for segmenting all of the clusters and
           %returning a label matrix
           
           data = bfopen([obj.saveLoc filesep obj.imageDataLoc]);
           
           im = loadStack(data,obj.imageDataInd);
           clear data;%Don't need to store this in memory anymore.
           
           [obj.labelm, ~, ~,~] = intensitySegmentLargeClusters(im,...
               obj.th,obj.minObjSize, obj.minObjVol, '10x volume');
           labelm = obj.labelm;
           save([obj.saveLoc filesep obj.labelSaveName], 'labelm');
       end
       
       function obj = segmentClusters_xds_dns(obj)
           %do we need this?
           %Run through the code for segmenting all of the clusters and
           %returning a label matrix
           
           data = bfopen([obj.saveLoc filesep obj.imageDataLoc]);
           
           im = loadStack(data,obj.imageDataInd);
           clear data;%Don't need this anymore.
           
           [obj.labelm, ~, ~,~] = intensitySegmentLargeClusters(im, obj.th,obj.minObjSize, obj.minObjVol, '10x volume xds dns');
           labelm = obj.labelm;
           save([obj.saveLoc filesep obj.labelSaveName], 'labelm');
       end
       
       function obj = segmentLCDClusters(obj)
           %Run through the analysis for identifying clusters in a LCD
           %background
           data = bfopen([obj.saveLoc filesep obj.imageDataLoc]);
           
           im = loadStack(data,obj.imageDataInd);
           clear data;%Don't need to store this in memory anymore
           
           [obj.labelm, ~, ~,~] = intensitySegmentLargeClusters(im, 100, 1000, 1000, '63x');
           labelm = obj.labelm;
           save([obj.saveLoc filesep obj.labelSaveName], 'labelm');
       end
       
       function obj = importData(obj)
           %Import the label matrix generated by segmentClusters
           inputvar = load([obj.saveLoc filesep obj.labelSaveName]);
           obj.labelm = inputvar.labelm;     
       end
       
       function obj = loadSegmentedData(obj)
           inputvar = load([obj.saveLoc filesep obj.saveName]);
           obj.cc = inputvar.obj.cc;
       end
       
       function obj = getRegionProps_2color(obj, im_1, im_2)
           %Collect image properties relevant to two color data sets. In
           %addition to calculating the area of each object, calculate the
           %the fractional occupancy of the aggregate by either color
           frac_thresh = 10;
           
           %Get properties of each identified object
           obj.cc = regionprops(obj.labelm,'Area', 'BoundingBox');
           obj = obj.filter;
          
           
           clear data;%Don't need to store this in memory anymore.
           
           %Get properties of each identified object
           obj.cc = regionprops(obj.labelm, im_1>frac_thresh, 'Area', 'BoundingBox', 'MeanIntensity');
           cc_temp = regionprops(obj.labelm, im_2>frac_thresh, 'Area', 'BoundingBox', 'MeanIntensity');
           %Calculate the fraction in each object
           for i=1:length(obj.cc)
              obj.cc(i).frac_1 = obj.cc(i).MeanIntensity;
              obj.cc(i).frac_2 = cc_temp(i).MeanIntensity;
           end
           obj.cc = rmfield(obj.cc, 'MeanIntensity');
       end
       function obj = getRegionProps(obj)
           %Get properties of each identified object
          obj.cc = regionprops(obj.labelm,'Area', 'BoundingBox');
          obj = obj.filter;    
       end
       
       function obj = loadRegionProps(obj)
          inputvar = load([obj.saveLoc filesep obj.saveName]);
          obj.cc = inputvar.obj.cc;
       end
       
       
       function obj = filter_LCD(obj)
           %Filter the LCD cluster data
           %As in filter function, pixel/micron conversion is specific to
           %the parameters that we use experimentally.
           
           for i=1:length(obj.cc)
               obj.cc(i).height = obj.cc(i).BoundingBox(end);
               obj.cc(i).volume = (1./10.7466)*(1./10.7466)*1*obj.cc(i).Area;
              obj.cc(i).ind = i;
           end
       end
       
       function obj = filter(obj)
           %Filter down the ist of clusters by removing clusters that don't
           %span many z-slices, or touch the boundary (maybe on the second)
           
           
           %MLJ: note need to replace hardcoded pixel size here with values
           %derived from imaging data-otherwise this is going to lead to
           %super-subtle bug!
           %if(length(obj.cc)==1 && isnan(obj.cc(1).volume))
          %    return; %Useful if this has already been allocated to be empty and we're running this code again for some reason 
          % end
               
           for i=1:length(obj.cc)
              obj.cc(i).height = obj.cc(i).BoundingBox(end);
              obj.cc(i).volume = 0.5679*0.5679*2*obj.cc(i).Area;
              obj.cc(i).ind = i;
           end
           %Put something there so that we can still go through the data
           %efficiently.
           
           if(isempty(obj.cc))
              obj.cc(1).volume = NaN;
              
           end
           %Only filter for actual objects, not place  holder values
           if(~isnan(obj.cc(1).volume))
               %Remove all clusters that don't span at least 3 z-slices
               ext = arrayfun(@(x)obj.cc(x).BoundingBox(end), 1:length(obj.cc));
               obj.cc(ext<2) = [];
           end
          %Remove all clusters with a volume smallers than 10^(2) ~ 5
          %microns cubed-two bacteria in size
          %  sz = arrayfun(@(x)obj.cc(x).volume, 1:length(obj.cc));
         % obj.cc(sz<5^3) = [];
           
    
       end
       
       function obj = constructHistogram(obj)
           %Construct a histogram of cluster volume
           if(isempty(obj.binsize)||isempty(obj.binlimits))
               [obj.histdata, obj.histedges] = histcounts([obj.cc.volume]);
           end
       end
       
       function obj = computeSummaryStatistics(obj)
           %Compute summary statistics for this distribution.   
           if(isempty(obj.cc))
               obj.median = NaN;
               obj.skewness = NaN;
           else
               obj.median = median([obj.cc.volume]);
               obj.skewness = skewness([obj.cc.volume]);
           end
       end
       
       function displayHistogram(obj)
           %Display histogram for this data
           %Make a dummy histogram and then fill it with the data generated
           %by the other function
           h = histogram(obj.histdata, obj.histedges);
           
          
           
       end
       
       
       function overlaySegmentation(obj)
           %Plot the result of the segmentation over the actual data.
           
           %Need to write this in two ways: One to briefly show the overlay
           %of the data, the other to allow the user to check whether the
           %segmentation worked succesfully or not.
           
           %Load image
           data = bfopen([obj.saveLoc filesep obj.imageDataLoc]);
           
           im = loadStack(data,obj.imageDataInd);
           clear data
           %Load mask
           inputvar = load(obj.labelSaveName);
           obj.labelm = inputvar.labelm;
           
           figure; 
           for i=1:size(obj.labelm,3)
               bw = obj.labelm(:,:,i)>0;
               bw = bwperim(bw);
               bw = imdilate(bw, strel( 'Disk', 1));
              
               thisim = im(:,:,i);
               imshow(thisim,[0 200]);
               hold on;
               green = cat(3, zeros(size(thisim)),ones(size(thisim)),zeros(size(thisim)));
               h = imshow(green);hold on;
               set(h, 'AlphaData', bw);
               
               title(num2str(i));
               pause
           end
           
           
       end
       
       
       function save(obj)
           save([obj.saveLoc filesep obj.saveName]);
       end
       
   end
   
end
    

classdef indiaInkQuantification
   properties
       thresh = 30;
       filename = '';
       im = []; %nx2 cell array, n = num of images to analyze, first index: fluorescent channel, second index: india ink
       imageName = {}; %Contains file names of each loaded image. Used to double check that we're loading everything in properly.
       prop = []; %Structured array containing information about a given data set.
       compiledData = {} %Structured array containing compiled data across replicates and summary statistics
   end
   
   methods
     
       
       
       function obj = loadStack(obj,filename, ind, varargin)
           %Load the image stack. The last field will always be the india
           %ink channel.
           obj.filename = filename;
           imAll = bfopen(filename);
               
           obj.im = cell(length(ind),2);
           for i=1:length(ind)
               obj.im{i,1} = imAll{ind(i),1}{1,1};
               obj.im{i,2} = imAll{ind(i),1}{2,1};
               
               obj.imageName{i} = imAll{ind(i),1}{1,2};
           end
       end
       
       function obj = calcProperties(obj)
           %Calculate the area and occupancy of each void given the defined threshold
           %Calculate mean and std
            for i=1:size(obj.im,1)
                segIm = obj.im{i,2}>obj.thresh;
                segIm = bwareaopen(segIm,1000);
                cc = regionprops(segIm, obj.im{i,1},'Area','PixelIdxList');
                if(isempty(cc))
                    obj.prop{i}.Area = [];
                    obj.prop{i}.Occupancy = [];
                    continue
                end
              for j=1:length(cc)
                 val = obj.im{i,1}(cc(j).PixelIdxList);
                 cc(j).Occupancy = sum(val(:)>0)/cc(j).Area;
              end
            obj.prop{i,1}.Area = [cc.Area];
            obj.prop{i,1}.Occupancy = [cc.Occupancy];
            end           
       end 
       
       function obj = loadMultiStack(obj, dirname, rootname, regnum, replnum, numcolor)
           %Load and analyze a series of images
           obj.filename = [dirname filesep rootname];
           
           %Temp. cludge to deal with sturcture in 8/21 data: brightfield
           %image in the second channel not third.
           clist = [1,2,3];
           for i = 1:regnum
               for j = 1:replnum
                   for c =1:numcolor
                       obj.im(:,:,c) = imread([obj.filename 'repl' num2str(j) '.tif'], 'Index', clist(c));
                   end
                   
                   %For each of the image sets, get properties of regions
                   %identified.ct
                   for c = 1:numcolor-1
                       segIm = obj.im(:,:,end)>obj.thresh;
                       segIm = bwareaopen(segIm,1000);
                       cc_temp = regionprops(segIm, obj.im(:,:,c),'Area','PixelIdxList'); 
                       if(isempty(cc_temp))
                           cc_temp.Area = [];
                           cc_temp.Occupancy = [];
                           obj.prop{i,j,c} = cc_temp;
                           continue
                       end
                       
                       for x=1:length(cc_temp)
                           imtemp = obj.im(:,:,c);
                           val = imtemp(cc_temp(x).PixelIdxList);
                           cc_temp(x).Occupancy = sum(val(:)>0)/cc_temp(x).Area;
                       end
                      
                       
                       obj.prop{i,j,c} = cc_temp;
                       
                   end
               end
           end
           
       end
   
       
       function obj = compileStatistics(obj)
          %Concatenate each of regions together and gather summary
          %statistics
          arraysize = size(obj.prop);
          if(length(arraysize)==2)
              %Update indez in case there is only one color channel
              arraysize(3) = 1;
          end
          
          for repl = 1:arraysize(2)
                  for c = 1:arraysize(3)
                      obj.compiledData.alldata{repl,c} = [];
                     for reg = 1:arraysize(1)
                         obj.compiledData.alldata{repl,c} = [obj.compiledData.alldata{repl,c} [obj.prop{repl,reg,c}.Occupancy]];
                     end
                  end
          end
          
          obj.compiledData.mean = cellfun(@(x)mean(x), obj.compiledData.alldata);
          
          
          
       end
   end

end
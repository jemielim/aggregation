%Analyze data from a time series of cluster data
% 
% <clusterAnalysisTimeSeries: Analysis of aggregate formation for a time series>
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
classdef clusterAnalysisTimeSeries
    
   properties
       saveLoc = '';
       timeList = [];
       numtp = [];
       genotype = '';
       alldata; %Location to save individual data sets
       
       %Summary statistics for time series
       medianval = [];
   end
   
   methods
       
       function obj = clusterAnalysisTimeSeries(input_saveLoc, input_timeList, input_genotype)
           obj.saveLoc = input_saveLoc;
           obj.timeList = input_timeList;
           obj.genotype = input_genotype;
           obj.numtp = length(obj.timeList);
           
           for t = 1:obj.numtp
               thisdata(t) = clusterAnalysis(obj.saveLoc, obj.genotype, obj.timeList(t));
           end
           
           obj.alldata = thisdata;
       end
       
       function obj = importAndFilterData(obj)
          for t = 1:obj.numtp
             obj.alldata(t) = obj.alldata(t).importData;
             obj.alldata(t) = obj.alldata(t).getRegionProps;
             obj.alldata(t).labelm = []; %Remove this so the files don't get too big.
          end
       end
       
       function obj = runPipeline(obj)
       %Run the segmentation analysis pipeline
           for i=1:length(obj.alldata)
             obj.alldata(i) = obj.alldata(i).runPipeline;
          end
       end
       
       function obj = reloadData(obj)
           %Load data from individual .mat files for each time points
           for t = 1:obj.numtp
               inputvar = load([obj.saveLoc filesep obj.alldata(t).saveName]);
               obj.alldata(t) = inputvar.obj;
           end
           
       end
       
       function obj = filter(obj)
           for i=1:obj.numtp
               obj.alldata(i) = obj.alldata(i).filter;
           end
       end
       
       function obj = constructHistogram(obj)
           %Construct a histogram of cluster volume
       end
       
       function obj = computeSummaryStatistics(obj)
           %Compute summary statistics for this distribution.
           for t = 1:obj.numtp
              obj.alldata(t) = obj.alldata(t).computeSummaryStatistics;
              obj.medianval(t) = obj.alldata(t).median;
           end
       end
       
       function save(obj)
           %Save object
           for i=1:obj.numtp
               save([obj.alldata(i).saveLoc filesep obj.alldata(i).saveName])
           end
       end
       
   end
   
end
    
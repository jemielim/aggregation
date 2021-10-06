%Analyze data from a set of time series data from different replicates and
%genotypes.
% <clusterReplicateAnalysis: Analysis of aggregate formation from groups of individual replicates>
%     Copyright (C) 2019, Matthew Jemielita
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
classdef clusterReplicateAnalysis
    
      
   properties
    allData = []; %file of data on all individual aggregates across all data sets
    normvol= 0.5679*0.5679*2*2048*2048*51;

    list = []; %List of all summary statistics from data set
    minVal = 10e-6;
    minNum = 20; %Minimum number of aggregates-if it's a small number, probably background noise.
    colorList = {'D47A', 'D47E', 'D47A vca0812', 'D47A vca0813', 'D47A vca0812_3'};
    colorVal = {[1 0 0], [1 0 1], [0 0 1], [0 0 0], [1 0.5 0]}
    genotypeList = {'D47A', 'D47E', 'D47A vca0812', 'D47A vca0813', 'D47A vca0812_3'};%List of all unique genotypes in this data set
    
    plotErrorBarHandle = []; %container for handles to errorbar plot so I can manipulate better
    offset = []; %Use for storing the offset from when D47A aggregates for each time series
    saveLoc = []; %Location to save this data file to
    
   end
   
   methods
       function obj = extract(obj,cts, repl)
       %Load in given cluster data and assign replicate number
       if(isempty(obj.allData))
           ind = 1;
       else
           ind = length(obj.allData)+1;
       end
       
       obj.allData(ind).repl = repl;
       obj.allData(ind).data = cts.alldata;
       obj.allData(ind).saveLoc = cts.saveLoc;
       obj.allData(ind).genotype = cts.genotype;
       end
       
       function obj = calcTotVol(obj)
           %Calculate the total volume
        
        for ind=1:length(obj.allData)
            for t= 1:length(obj.allData(ind).data)
                if(isempty(obj.allData(ind).data(t).cc))
                   obj.allData(ind).totvol(t) = obj.minVal; %Minimum set to make plotting easier
                elseif(isnan([obj.allData(ind).data(t).cc.volume]))
                    obj.allData(ind).totvol(t) = obj.minVal; %Minimum set to make plotting easier
                else
                    
                    %Minimum number of objects in 
                    if(length([obj.allData(ind).data(t).cc.volume])>obj.minNum)
                        obj.allData(ind).totvol(t) = sum([obj.allData(ind).data(t).cc.volume])/obj.normvol;
                    else
                        obj.allData(ind).totvol(t) = obj.minVal;
                    end
                    
                    %if equal to zero, set to minimum to enable plotting
                    if(obj.allData(ind).totvol(t) ==0)
                        obj.allData(ind).totvol(t) = obj.minVal;
                    end
                end
            end
        end
        
       end 
       function obj = calcAverageSize(obj)
           %Calculate average aggregate size (vol*vol/totvol: normalized by
           %average aggregates size a cell is in)
           for ind=1:length(obj.allData)
               for t= 1:length(obj.allData(ind).data)
                   if(isempty(obj.allData(ind).data(t).cc))
                       obj.allData(ind).aveSize(t) = 0;
                   elseif(isnan([obj.allData(ind).data(t).cc.volume]))
                       obj.allData(ind).aveSize(t) = 0;
                   else
                       
                       %Minimum number of objects in
                       if(length([obj.allData(ind).data(t).cc.volume])>obj.minNum)
                           obj.allData(ind).aveSize(t) = sum([obj.allData(ind).data(t).cc.volume].*[obj.allData(ind).data(t).cc.volume])./sum([obj.allData(ind).data(t).cc.volume]);
                       else
                           obj.allData(ind).aveSize(t) = 0;
                       end
                       
                   end
               end
           end
        
       end
           
       function obj = getMasterList(obj)
           %Assemble all master data into one list so that it's easier to
           %assemble summary statistics and plots, etc.
           
           for i=1: length(obj.allData)
               for j=1:length(obj.allData(i).data)
               obj.list.totvol(i,j) = [obj.allData(i).totvol(j)];
               obj.list.avesize(i,j) = [obj.allData(i).aveSize(j)];
               obj.list.genotype{i,j} = obj.allData(i).data(j).genotype;
               obj.list.time(i,j) = obj.allData(i).data(j).time;
               
               
               repl = obj.allData(i).repl;
               ofs = obj.offset(repl);
               thisTime = [obj.allData(i).time]-obj.allData(i).time(ofs);
               obj.list.offsetTime(i,j) = thisTime(j);
               
               end
           end
           
           
       end
       function obj = getReplicateError(obj)
           %Find the mean and std dev for each time point (offset) across
           %replicates of the same genotype
           
          for i=1:length(obj.genotypeList)
              %Get indices corresponding to each genotype
              ind = cellfun(@(x)strcmp(x, obj.genotypeList{i}), obj.list.genotype, 'UniformOutput', false);
              ind = cell2mat(ind);
              ind = find(ind==1);
              indGenotype{i} = ind;
              
              offsetTime = obj.list.offsetTime(ind);
              totvol = obj.list.totvol(ind);
              uniqt = unique(offsetTime);
              
              for j=1:length(uniqt)
                 indT = find(offsetTime == uniqt(j));
                 obj.list.mean(i,j)= mean([totvol(indT)]);
                 obj.list.stdDev(i,j)= std([totvol(indT)]);
                 obj.list.timeReplicate(i,j) = uniqt(j);
              end

          end
            
       end
       
       function obj = getReplicateErrorAveSize(obj)
           %Find the mean and std dev for average aggregate size at each
           %time point
           for i=1:length(obj.genotypeList)
              %Get indices corresponding to each genotype
              ind = cellfun(@(x)strcmp(x, obj.genotypeList{i}), obj.list.genotype, 'UniformOutput', false);
              ind = cell2mat(ind);
              ind = find(ind==1);
              indGenotype{i} = ind;
              
              offsetTime = obj.list.offsetTime(ind);
              avesize = obj.list.avesize(ind);
              uniqt = unique(offsetTime);
              
              for j=1:length(uniqt)
                 indT = find(offsetTime == uniqt(j));
                 obj.list.sizemean(i,j)= mean([avesize(indT)]);
                 obj.list.sizestdDev(i,j)= std([avesize(indT)]);
                 obj.list.timeReplicate(i,j) = uniqt(j);
                 
                 %set to nan if zero to make plotting easier
                 if(mean([avesize(indT)])==0)
                      obj.list.sizemean(i,j)= NaN;
                      obj.list.sizestdDev(i,j) = NaN;
                 end
              end

          end
       end
       
       function obj = getTime(obj)
          %Add time list to each data set
          for i=1:length(obj.allData)
             t = [];
             for j=1:length(obj.allData(i).data)
                t = [t obj.allData(i).data(j).time];
             end
             obj.allData(i).time = t;
          end
          
       end
       
       function obj = calcOffset(obj)
          %Calculate the offset from the time that D47A aggregates for each
          %time series
          replList = [];
          genList = [];
          
          for i=1:length(obj.allData)
              ind = obj.allData(i).totvol>0.001; ind = find(ind,1);
              obj.allData(i).offset = obj.allData(i).time(ind);
              if(string(obj.allData(i).data(1).genotype) == 'D47A')
                  obj.offset(obj.allData(i).repl) = ind;
              end
              
          end 
          
       end
       
       function obj = runPipeline(obj)
           %Run all modules needed for analysis
          obj = obj.calcTotVol;
          obj = obj.calcAverageSize;
          obj = obj.getTime;
          obj = obj.calcOffset;
          
          obj = obj.getMasterList;
          
          obj = obj.getReplicateError;
          obj = obj.getReplicateErrorAveSize;
       end
       
       function obj = plotTotVol(obj)
           
           figure; hold on
           for i=1:length(obj.allData)
               
               cl = cellfun(@(x)isequal(obj.allData(i).data(1).genotype, x), obj.colorList);
               
               ind = find(cl==1);
               
               plot(obj.allData(i).totvol, 'Color', obj.colorVal{ind});
           end
           set(gca, 'YScale', 'log')
           set(gca, 'YLim', [0.1*obj.minVal, 1]);
       end
    
       function obj = plotTotVolErrorBar(obj) 
           %Plot all data sets using mean and std. dev.
           
           figure; hold on
           for i=1:size(obj.list.mean,1)
               x = obj.list.timeReplicate(i,:);
               y = obj.list.mean(i,:);
               sd = obj.list.stdDev(i,:);
               h(i) = plot(x,y,'Color', obj.colorVal{i}, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 24);
               p(i) = errorbar(x,y,sd, 'Color', obj.colorVal{i});
               set(p(i), 'LineWidth', 1.5);
           end
           set(gca, 'YScale', 'log')
           set(gca, 'YLim', [0.1*obj.minVal, 1]);
           set(gca,'XLim', [x(1), x(end)+1]);
           
           
           l(3) = legend(h, obj.colorList, 'Location', 'Northwest');

           l(1) = xlabel('Time from \DeltavpsL HCD-locked aggregation onset (h)');
           l(2) = ylabel('Aggregate volume fraction');
           
           set(gca, 'FontSize',14)
           set(gca, 'LineWidth',1.5);
           
           %Pass out handles so I can modify this figure readily
           obj.plotErrorBarHandle.h = h;
           obj.plotErrorBarHandle.p = p;
           obj.plotErrorBarHandle.l = l;
           obj.plotErrorBarHandle.a = gca;
       end
      
       function obj = plotTotVolGenotype(obj, genotype)
           %Plot only data for one genotype
           figure; hold on
           n = 1;
           for i=1:length(obj.allData)
               if(string(obj.allData(i).data(1).genotype)==genotype)
                   h(i) = plot(obj.allData(i).totvol);
                   leg{n} = obj.allData(i).genotype;n=n+1;
               end
           end
           legend(leg);
           set(gca, 'YScale', 'log')
           set(gca, 'YLim', [0.1*obj.minVal, 1]);

           title(genotype)
           
       end
       
       function obj = plotTotVolOffset(obj)
           %Plot total volume with offset with respect to when D47A starts
           %to aggregate for giving time series
           figure; hold on
           n = 1;
           
           for i=1:length(obj.allData)
               
               %Calculate time offset from onset of aggregation in D47A
               %strain
               repl = obj.allData(i).repl;
               ofs = obj.offset(repl);
               thisTime = [obj.allData(i).time]-obj.allData(i).time(ofs);
               
               %Get color map    
               cl = cellfun(@(x)isequal(obj.allData(i).data(1).genotype, x), obj.colorList);
               ind = find(cl==1);
               
               
               h(i) = plot(thisTime,obj.allData(i).totvol,'Color', obj.colorVal{ind}, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 24);
               l(i) = ind;
           end
           
           %Find first unique index for each color-otherwise you get a lot
           %of overlap in colors used
           ind = arrayfun(@(x)find(l==x,1, 'first'), 1:length(obj.colorList));
           legend(h(ind), obj.colorList(l(ind)), 'Location', 'Northwest');
           
           set(gca, 'YScale', 'log')
           set(gca, 'YLim', [0.1*obj.minVal, 1]);
           xlabel('Time from HCD-locked aggregation onset (h)')
           ylabel('Aggregate volume fraction');
           
           set(gcf, 'Position', 1.0e+03*[4.1050 0.7743 0.7247 0.5447]);
           set(gca, 'Position', [0.13 0.11 0.7750 0.8150]);
           set(gca, 'FontSize',14)
           set(gca, 'LineWidth',1.5);
       end
       
        function obj = plotAveSizeErrorBar(obj) 
       %Plot all data sets using mean and std. dev.
       
        figure; hold on
           for i=1:size(obj.list.mean,1)
               x = obj.list.timeReplicate(i,:);
               y = obj.list.sizemean(i,:);
               sd = obj.list.sizestdDev(i,:);
               h(i) = plot(x,y,'Color', obj.colorVal{i}, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 24);
               p = errorbar(x,y,sd, 'Color', obj.colorVal{i});
               set(p, 'LineWidth', 1.5);
           end
           set(gca, 'YScale', 'log')
     
           set(gca,'XLim', [x(1), x(end)+1]);
           
           
           legend(h, obj.colorList, 'Location', 'Northwest');

           xlabel('Time from HCD-locked aggregation onset (h)')
           ylabel('Aggregate volume (\mum)^3');
           
           set(gcf, 'Position', 1.0e+03*[4.1050 0.7743 0.7247 0.5447]);
           set(gca, 'Position', [0.13 0.11 0.7750 0.8150]);
           set(gca, 'FontSize',14)
           set(gca, 'LineWidth',1.5);
           
          
        end
       
       
       function obj = save(obj)
           %Save file
           save([obj.saveLoc filesep 'clusterReplicate.mat'], 'obj');
       end
       
       
   end

   
end
      
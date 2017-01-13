files2run = {'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10262016_alldays.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10282016_alldays.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10312016_alldays.mat'};
         
 [PO,PDOD1,PDOD2] = deal(cell(length(files2run),1));
 
 for fi = 1:length(PO) 
     load(files2run{fi});
     Stay_switch_check;
     PO{fi} = percover;
     clearvars -except PO PDOD1 PDOD2 files2run
 end
 
 for i = 1:length(PO)
     PDOD1{i} = [PO{i}.pd{1}; PO{i}.od{1}];
     PDOD2{i} = [PO{i}.pd{2}; PO{i}.od{2}];
 end
 
 PDODall{1} = cell2mat(PDOD1);
 PDODall{2} = cell2mat(PDOD2);
 %%
 PDall{1} = cell2mat(cellfun(@(x) sum(x.pd{1},2),PO,'Uni',0));
 ODall{1} = cell2mat(cellfun(@(x) sum(x.od{1},2),PO,'Uni',0));
 
 PDall{2} = cell2mat(cellfun(@(x) sum(x.pd{2},2),PO,'Uni',0));
 ODall{2} = cell2mat(cellfun(@(x) sum(x.od{2},2),PO,'Uni',0));
 
 PDall{1}(isnan(PDall{1})) = [];
 PDall{2}(isnan(PDall{2})) = [];
 ODall{1}(isnan(ODall{1})) = [];
 ODall{2}(isnan(ODall{2})) = [];
%  
%  figure; hold on; 
%  linehist(linspace(0,100,50),PDODall{1}(:,1));
%  linehist(linspace(0,100,50),PDODall{2}(:,1));
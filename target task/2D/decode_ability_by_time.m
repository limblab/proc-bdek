%% Set parameters for neural binning
t12T = [-.2 0.7]; % in seconds
t12G = [-.2 0.4]; % in secconds
bsize = 50; % in milliseconds

Tis = 1:round(diff(t12T)*1000/bsize);
Gis = (1:round(diff(t12G)*1000/bsize))+Tis(end);
bl_i = 1:sum(t12T(1):bsize/1000:t12T(end) < 0);

tt = alldays(1).tt;
units = alldays(1).PMd_units;
%% Loop through units and create rasters
[raster, rasterB] = deal(cell(length(units),1));
for i = 1:length(units)
    clc; fprintf('%d/%d\n',i,length(units));
    rasterT = raster_get(units{i}(2:end),tt,...
        t12T,'target'); % Create target aligned rasters
    rasterG = raster_get(units{i}(2:end),tt,...
        t12G,12); % Create movement aligned rasters
    raster{i} = [rasterT rasterG]; % Concatenate
    rasterB{i} = bin_array(raster{i},size(raster{i},1),...
                          (diff(t12T)+diff(t12G))*1000/bsize,'sum');
end

%% Categorize reach directions
C = mod(round(tt(:,2)/(pi/4)), 8)+1;


%%
av8 = zeros(8,size(rasterB{1},2));
maxmod = zeros(length(rasterB),size(rasterB{1},2));

avcount = round(size(tt,1)/8);
randfunc = @(x) abs(diff(nanmean(reshape(x(randperm(length(x),avcount*2)),[],2))));
bdist = cell(length(rasterB),1);
for n = 1:length(rasterB) % n = neuron
    
    blines = reshape(rasterB{n}(:,bl_i),[],1);
    blinemean = nanmean(blines);
    
    for d = 1:8 % d = target direction 
        av8(d,:) = nanmean(rasterB{n}(C==d,:));
    end
     
    maxmod(n,:) = max(av8)-min(av8);
    
    for t = 1:size(av8,2)

       [av8L, av8H] = boot_bounds(1000,@nanmean,rasterB{n}(C==d,t),2.5,97.5); 
       if (blinemean > av8L) && (blinemean < av8H)
           sigmat(n,t
    end

%     [~,~,~,bdist{n}] = boot_bounds(1000,randfunc,blines,2.5,97.5);
    clc; fprintf('%d/%d\n',n,length(rasterB));
end

%%
p_of_C = zeros(length(rasterB),size(rasterB{1},2));
for n = 1:length(rasterB)

    [val,ord] = sortrows(maxmod(n,:)');
end
    

%%
% %% Build decoder
% AVS = cell(length(rasterB),1);
% for i = 1:length(rasterB) % Loop through units
%     r = rasterB{i};
%     clc; fprintf('%d/%d\n',i,length(rasterB));
%     for t = 1:size(rasterB{1},2) % Loop through time points
%         for j = 1:size(unique(C),1) % loop through targets
% 
%             AVS{i}(j,t) = nanmean(r(C==j));
%         end
%     end
% end

%%
% decode = cell(length(rasterB),1);
% for i = 1:length(rasterB) % i = unit
%     
%     r = rasterB{i};
%     clc; fprintf('%d/%d\n',i,length(rasterB));
%     for j = 1:size(rasterB{1},1) % j = trial
%        
%         targ = C(j);            
%         count = r(j,:);
%         av_same = nanmean(r(C==targ & ~ismember((1:size(rasterB{1},1))',j),:));
%         av_others = AVS{i}(~ismember(1:8,targ),:);
%         
%         ptarg_same = poisspdf(count,av_same);
%         ptarg_others = poisspdf(count,av_others);
%         for t = 1:size(rasterB{1},2) % t = time
% 
% 
% 
%             decode{i}(j,t) = ptarg_same > max(ptarg_others);
% 
%         end
%     end
% end
%             
%         
        
        
        
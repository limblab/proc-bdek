function[cents,dirbinrast,dirbinlow,dirbinhigh,ind,av_spiketime,full_rasts,full_rasts_inds] = co_tuning(spikes,angles,comp_tt,t1,t2,alignment,neuron,bootsamps,varargin)

direcs = angles;
direcs(direcs < 0) = direcs(direcs < 0) + 2*pi;
cents = 0:pi/4:(2*pi-pi/4); 

train = spikes{neuron};
like_conds = unique(comp_tt(:,3));

if strcmp(alignment,'target')
    aligntype = 5;
elseif strcmp(alignment,'go')
    aligntype = 6;
elseif strcmp(alignment,'end')
    aligntype = 7;
elseif isnumeric(alignment)
    aligntype = alignment;
else
    fprintf('alignment not supported: using target alignment...\n');
end

raster = zeros(length(comp_tt),t2-t1);
avtime = zeros(length(comp_tt),1);
for i = 1:length(comp_tt)
    
    timestart = comp_tt(i,aligntype)+t1./1000;
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=(t2-t1)) = [];
    
    avtime(i) = mean(aligned_ts);
    
    raster(i,aligned_ts) = 1;
    
end

rastfull = raster(:,1:(t2-t1));
rastbin = sum(rastfull,2);
rastbin = rastbin./((t2-t1)/1000);

ind = cell(8,1);
[dirbinrast, dirbinlow, dirbinhigh, av_spiketime] = deal(zeros(8,1));
full_rasts = cell(8,length(like_conds));
full_rasts_inds = cell(8,length(like_conds));

if nargin < 8
    bootsamps = 0;
elseif isempty(bootsamps)
    bootsamps = 0;
end

for i = 1:8
    dists_from_cent = angle_diff(direcs,repmat(cents(i),length(direcs),1));
    ind{i} = find(abs(dists_from_cent) < pi/4);
    dirbinrast(i) = nanmean(rastbin(ind{i}));
    av_spiketime(i) = nanmean(avtime(ind{i}));
    
    for j = 1:length(like_conds)
        
        likeinds = find(comp_tt(:,3)==like_conds(j));
        
        dir_like_inds = ind{i}(ismember(ind{i},likeinds));
        
        if ~isempty(dir_like_inds)
            full_rasts{i,j} = rastfull(dir_like_inds,:);
            full_rasts_inds{i,j} = dir_like_inds;
        end
    end
    
    
    if bootsamps > 1
        if isnan(dirbinrast(i)) || isempty(dirbinrast(i))
            dirbinlow(i) = NaN;
            dirbinhigh(i) = NaN;
        elseif length(ind{i}) == 1
            dirbinlow(i) = rastbin(ind{i});
            dirbinhigh(i) = rastbin(ind{i});
        else
            [dirbinlow(i), dirbinhigh(i)] = boot_bounds(bootsamps,@nanmean,rastbin(ind{i}),2.5,97.5);
        end
    else
        dirbinlow(i) = dirbinrast(i);
        dirbinhigh(i) = dirbinrast(i);
    end
end


end

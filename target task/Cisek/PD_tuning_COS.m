function[PD,PD_bounds,PD_dists,MAG_dists,curve_params] = PD_tuning_COS(spikes,angles,comp_tt,t1,t2,alignment,neuron,speeds,bootstrapnum,bounds,varargin)

direcs = angles;
direcs(direcs < 0) = direcs(direcs < 0) + 2*pi;

train = spikes{neuron};

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
for i = 1:length(comp_tt)
    
    timestart = comp_tt(i,aligntype)+t1./1000;
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=(t2-t1)) = [];

    raster(i,aligned_ts) = 1;
    
end
rastfull = raster(:,1:(t2-t1));
rastbin = sum(rastfull,2)./((t2-t1)/1000);

if ~isempty(speeds)
    input_mat = [direcs rastbin speeds];
else
    input_mat = [direcs rastbin];
end

[~,curve_params] = PD_helper_func_COS(input_mat);
if nargin > 8
    [~,~,~,rand_boot] = boot_bounds(bootstrapnum, @PD_MAG_helper_func_COS, input_mat,bounds(1),bounds(2));
    PD = circ_mean(rand_boot(:,1));
    delta_from_av = circ_dist(rand_boot(:,1),PD);
    PD_bounds = prctile(delta_from_av,bounds);  
    PD_dists = circ_dist(rand_boot(:,1),PD);
    MAG_dists = rand_boot(:,2);%prctile(sort_boot(:,2),bounds);
else
    PD_bounds = [nan nan];
    MAG_dists = [nan];
    PD_dists = [nan];
end


end

function[cents,dirbinrast] = co_tuning_boot(spikes,angles,comp_tt,t1,t2,alignment,neuron,bootsamps)

direcs = angles;
direcs(direcs < 0) = direcs(direcs < 0) + 2*pi;
cents = 0:pi/4:(2*pi-pi/4); 

train = spikes{neuron};

if strcmp(alignment,'target')
    aligntype = 5;
elseif strcmp(alignment,'go')
    aligntype = 6;
else
    fprintf('alignment not supported: using target alignment\n');
end

raster = zeros(length(comp_tt),t2-t1);

for i = 1:length(comp_tt)
    
    timestart = comp_tt(i,aligntype)+t1./1000;
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=(t2-t1)) = [];
    
    raster(i,aligned_ts) = 1;
    
end

rastbin = sum(raster(:,1:(t2-t1)),2);
rastbin = rastbin./((t2-t1)/1000);

ind = cell(8,1);
checkrast = zeros(8,1);
dirbinrast = zeros(bootsamps,8);

for i = 1:8
    dists_from_cent = angle_diff(direcs,repmat(cents(i),length(direcs),1));
    ind{i} = find(abs(dists_from_cent) < pi/4);
    checkrast(i) = nanmean(rastbin(ind{i}));
    
    if isempty(checkrast(i))
        dirbinrast(:,i) = NaN;
    else
        dirbinrast(:,i) = bootstrp(bootsamps,@nanmean,rastbin(ind{i}));
    end

end


end

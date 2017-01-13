function[thetas,curve,p,curveL,curveH] = co_tuning_VM(spikes,angles,comp_tt,t1,t2,alignment,neuron,bin,bootstrapnum,varargin)

direcs = angles;
direcs(direcs < 0) = direcs(direcs < 0) + 2*pi;
cents = 0:pi/4:(2*pi-pi/4); 

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
rastbin = sum(rastfull,2);
rastbin = rastbin./((t2-t1)/1000);

ind = cell(8,1);
[dirbinrast] = deal(zeros(8,1));
thetas = 0.01:0.01:2*pi;

if bin ==1
    for i = 1:8
        dists_from_cent = angle_diff(direcs,repmat(cents(i),length(direcs),1));
        ind{i} = find(abs(dists_from_cent) < pi/4);
        dirbinrast(i) = nanmean(rastbin(ind{i}));
    end
    
    [p,~,curve] = VM_fit3(cents,dirbinrast,thetas);
    if nargin > 8
        [curveL,curveH] = boot_bounds(bootstrapnum, @VM_helper_func, [cents',dirbinrast],2.5,97.5);
    else
        [curveL,curveH] = deal(NaN(1,1));
    end

else
    [p,~,curve] = VM_fit3(direcs,rastbin,thetas);
    if nargin > 8
        [curveL,curveH] = boot_bounds(bootstrapnum, @VM_helper_func, [direcs,rastbin],2.5,97.5);
    else
        [curveL,curveH] = deal(NaN(1,1));
    end
end

%vm_func = @(xs,p) p(1) + p(2)*exp(p(3)*cos(xs-p(4)));

%curve = vm_func(thetas,p);


end

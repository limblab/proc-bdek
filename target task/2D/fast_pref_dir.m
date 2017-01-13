function[prefd,cents,dirbinrast,dirbinstd,rasters] = fast_pref_dir(spikes,angles,comp_tt,nbins,t1,t2,alignment,neuron,varargin)

do_plot = 0;

direcs = angles;
direcs(direcs < 0) = direcs(direcs < 0) + 2*pi;
edges = linspace(0,2*pi,nbins+1);
cents = .5*(edges(2:end)+edges(1:end-1));

if nargin > 7
    
    train = spikes{neuron};

    [rast,rastbin] = bin_rast(comp_tt,train,t1,t2,1,alignment);
    rastbin = rastbin./((t2-t1)/1000);

    ind = cell(nbins,1);
    dirbinrast = zeros(nbins,1);
    dirbinstd = zeros(nbins,1);
    rasters = nan(1000,nbins);
    for i = 1:nbins
        ind{i} = find(direcs > edges(i) & direcs < edges(i+1));
        dirbinrast(i) = mean(rastbin(ind{i}));
        dirbinstd(i) = std(rastbin(ind{i}))./sqrt(length(rastbin(ind{i})));
        rasters(1:length(ind{i}),i) = rastbin(ind{i}).*((t2-t1)/1000);
    end

    if do_plot==1
        figure; errorbar(cents,dirbinrast,dirbinstd,'b.-');
    end
    prefd = [cents(dirbinrast==max(dirbinrast)) min(dirbinrast) max(dirbinrast)];



else
    
    prefd = zeros(length(spikes),4);
    rastbin = cell(length(spikes),1);
    figure; hold on;
    plot(cos(0:0.1:2*pi),sin(0:0.1:2*pi),'k');
    for neu = 1:length(spikes)
        fprintf('neuron: %d/%d\n',neu,length(spikes));
        
        train = spikes{neu};
        
        [rast,rastbin{neu}] = bin_rast(comp_tt,train,t1,t2,1,alignment);

        ind = cell(nbins,1);
        dirbinrast = zeros(nbins,1);
        for i = 1:nbins
            ind{i} = find(direcs > edges(i) & direcs < edges(i+1));
            dirbinrast(i) = mean(rastbin{neu}(ind{i}));
        end
        
        peakloc = find(dirbinrast==max(dirbinrast));
        mod_depth = (max(dirbinrast)-min(dirbinrast))./max(dirbinrast);
        
        prefd(neu,:) = [cents(peakloc(1)),min(dirbinrast),max(dirbinrast),mod_depth];
        
        plot([0 mod_depth.*cos(prefd(neu,1))],[0 mod_depth.*sin(prefd(neu,1))],'b');
        plot(mod_depth.*cos(prefd(neu,1)),mod_depth.*sin(prefd(neu,1)),'g.');
        
        drawnow; clc;
    end

end

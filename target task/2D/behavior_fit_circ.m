function [PLR,PLRb,kratiofit,priors,fitlines,resids] = behavior_fit_circ(targ_cent_reach_label, set_krats, set_priors, varargin)

targ_cent_reach_label(isnan(sum(targ_cent_reach_label,2)),:) = [];

targets = targ_cent_reach_label(:,1);
centroids = targ_cent_reach_label(:,2);
reach = targ_cent_reach_label(:,3);
labels = targ_cent_reach_label(:,4);

% mufunc = @(mu1,mu2,kratio) mu1 + atan(sin(mu2-mu1)./(kratio + cos(mu2-mu1)));
mufunc = @(mu1,mu2,kratio) mu1 + atan2(sin(mu2-mu1),(kratio + cos(mu2-mu1)));
mufitfunc = @(mu1,mu2,reaches,kratio) sum(circ_dist(mu1 + atan(sin(mu2-mu1)./(kratio + cos(mu2-mu1))),reaches).^2);
k2func = @(mu1,mu2,k3,kratio) sqrt(k3^2./(kratio^2 + 2*kratio*cos(mu2-mu1) + 1));

ls = flipud(unique(labels));
[kratfit,estreaches,PLR,PLRb,is,kratiofit,fitlines,resids] = deal(cell(length(ls),1));
priors = zeros(1,length(ls));
for i =  1:length(ls)
    is{i} = find(labels==ls(i));
    
    if nargin > 2
        mu1 = repmat(set_priors(i),length(is{i}),1);
    else
        mu1 = repmat(circ_mean(targets(is{i})),length(is{i}),1);
    end
    priors(i) = unique(mu1);
    
    mu2 = centroids(is{i});
    reaches = reach(is{i});
    
    if nargin > 1
        kratfit{i} = set_krats{i};
    else
        kratfit{i} = fminbnd(@(x) mufitfunc(mu1,mu2,reaches,x),eps,1000);
    end
    repkrat = repmat(kratfit{i},length(is{i}),1);

    estreaches{i} = mufunc(mu1,mu2,kratfit{i});
    
    k3 = circ_kappa(circ_dist(estreaches{i},reaches));
    k2 = mean(k2func(mu1,mu2,k3,kratfit{i}));
    k1 = kratfit{i}*k2;

    PLR{i} = [k1 k2 k3]; PLR{i}(PLR{i}<0.5)=.5;
    
    kratiofit{i} = PLR{i}(1)./PLR{i}(2);
    
    fitlines{i}(1,:) = min(centroids):0.01:max(centroids);
    fitlines{i}(2,:) = mufunc(priors(i),(fitlines{i}(1,:))',kratfit{i});
    
    [~,~,~,PLRb{i}] = boot_bounds(1000,@behavior_fit_circ_helper,[mu1,mu2,reaches,repkrat],2.5,97.5);
    
    resids{i} = circ_dist(estreaches{i},reaches);
end

end
                      
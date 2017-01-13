function [PLR] = behavior_fit_func(centroids_reaches_labels)

centroids = centroids_reaches_labels(:,1);
reaches = centroids_reaches_labels(:,2);
labels = centroids_reaches_labels(:,3);

ls = flipud(unique(labels));
% [kapres,dat_slope,KL,KP,KR,sigR,sigP,sigL,truesigP,truesigL,truesigR] = deal(zeros(length(ls),1));
[doi,rs,PLR] = deal(cell(length(ls),1));

bootfunc = @(res_slope) [1/sqrt((1-unique(res_slope(:,2)))*circ_kappa(res_slope(:,1))), ...
                         1/sqrt(unique(res_slope(:,2))*circ_kappa(res_slope(:,1))), ...
                         1/sqrt(circ_kappa(res_slope(:,1)))]; 

for i =  1:length(ls)
    is = find(labels==ls(i));
    
    [doi{i},rs{i}] = circ_polyfit(centroids(is),reaches(is));
    
    resslope = [rs{i},repmat(doi{i}(1),size(rs{i},1),1)];
    
    [~,~,PLR{i}] = boot_bounds(1000,bootfunc,resslope,2.5,97.5);
    
    PLR{i}(imag(PLR{i})~=0) = inf;
end



end
                      
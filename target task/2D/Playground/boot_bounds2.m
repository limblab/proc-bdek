function[lower_bound, upper_bound, sort_boot, rand_boot] = boot_bounds2(nboot,bootfun,data,lower_perc,upper_perc,dim,varargin)

if upper_perc < 1
    lower_perc = lower_perc.*100;
    upper_perc = upper_perc.*100;
end


booted = bootstrp(nboot,bootfun,data);%,'Options',options);

lower_bound = zeros(size(booted,2),1);
upper_bound = zeros(size(booted,2),1);

[sort_boot,rand_boot] = deal(zeros(nboot,size(booted,2)));
for i = 1:size(booted,2)
    sort_boot(:,i) = sortrows(booted(:,i));
    rand_boot(:,i) = booted(:,i);
    
    lower_bound(i) = sort_boot(round(lower_perc*size(sort_boot(:,i),1)),i);
    upper_bound(i) = sort_boot(round(upper_perc*size(sort_boot(:,i),1)),i);
end

end


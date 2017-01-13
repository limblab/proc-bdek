%'PMd SN  6251-000987.cmp'
%'M1 SN  6250-000989.cmp'
map = ArrayMap('PMd SN  6251-000987.cmp');
unit_activity = PMd_units;

unit_metrics = PMd_cdisc;
unit_categories = PMd_cdisc_ind;

unit_ids = zeros(length(unit_activity),1);
for i = 1:length(unit_activity)
    unit_ids(i) = unit_activity{i}(1);
end
unit_ids(end-1:end) = [];

metric_map = zeros(size(map,1),size(map,2),max(round(10*(unit_ids - floor(unit_ids)))));
categ_map = zeros(size(map,1),size(map,2),max(round(10*(unit_ids - floor(unit_ids)))));
flipped_map = metric_map;
flipped_cat = categ_map;
for i = 1:length(unit_ids)
    
    channel = floor(unit_ids(i));
    unitnum = round(10*(unit_ids(i)-floor(unit_ids(i))));
    
    [locx,locy] = find(map==channel);
    
    metric_map(locx,locy,unitnum) = unit_metrics(i);
    categ_map(locx,locy,unitnum) = unit_categories(i);
end

for i = 1:size(metric_map,3)
    flipped_map(:,:,i) = flipud(metric_map(:,:,i));
    flipped_cat(:,:,i) = flipud(categ_map(:,:,i));
end
    
plot_colors = {'c','b','k'};
figure; hold on;
for k = 1:size(metric_map,3)
    for j = 1:size(metric_map,2)
        for i = 1:size(metric_map,1)
            if flipped_map(i,j,k) ~= 0
                if flipped_map(i,j,k) > pi/4
                    plot(j,i,'o','Color',plot_colors{flipped_cat(i,j,k)},'MarkerSize',200*(flipped_map(i,j,k) - pi/4),'LineWidth',4);
                else
                    plot(j,i,'r.','MarkerSize',20,'LineWidth',4);
                end
            end
        end
    end
end
axis square
axis off
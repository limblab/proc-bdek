area = 'PMd';

baseline = zeros(length(eval(sprintf('alldays(1).%s_units',area))),1);
for i = 1:length(eval(sprintf('alldays(1).%s_units',area)))
    
    clc; fprintf('%d/%d\n',i,length(eval(sprintf('alldays(1).%s_units',area))));
    r = raster_plot(eval(sprintf('alldays(1).%s_units{%d}',area,i)),alldays(1).tt,[0 0.8],'center',0,'none');
    
    r_av = sum(sum(bin_array(r{1},size(r{1},1),8,'sum')))./(8*size(r{1},1))./(0.1);
    
    baseline(i,:) = r_av;
end
    
    
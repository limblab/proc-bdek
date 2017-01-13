
units = load('C:\Users\limblab\Desktop\Mihili_bdf\Unit_info\Mihili_PMd_Units_07192013.mat');
%%
unitcell = struct2cell(units);
waves = vertcat(cell2mat(cellfun(@(x) x(4:51),unitcell,'UniformOutput',0)));

norm_waves = (waves-repmat(min(waves,[],2),1,size(waves,2)))./repmat(max(waves,[],2)-min(waves,[],2),1,size(waves,2));

trough2peak = zeros(size(waves,1),1);
for i = 1:size(waves,1)
    trough2peak(i) = find(norm_waves(i,:)==1,1,'first')-find(norm_waves(i,:)==0,1,'last');
end

n_w = kmeans(trough2peak,2);

narrow_inds = find(n_w == 1);
broad_inds = find(n_w == 2);
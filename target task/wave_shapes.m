function[waves,wave_bounds] = wave_shapes(bdf)

ids = vertcat(bdf.units.id);

sorted_inds = find(ids(:,2) ~= 0 & ids(:,2) ~= 255);

[waves,wave_bounds] = deal(cell(length(sorted_inds),1));
for i = 1:length(sorted_inds)
        
    wave = bdf.units(sorted_inds(i)).wave(1,:);
    stdwave = bdf.units(sorted_inds(i)).wave(2,:);
    
    waves{i} = wave;
    
    wave_bounds{i}(1,:) = wave + stdwave;
    wave_bounds{i}(2,:) = wave - stdwave;

end

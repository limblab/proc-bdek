modaye = [MO DA YE];
orig_place = cd;

%% Go to directory
if strcmp(computer,'work')
    cd(sprintf('C:\\Users\\limblab\\Desktop\\%s_bdf\\Unit_info\\',monkey));
else 
    fprintf('Manually Load File\n');
end

%% If it exists, load the file
variables = what;
if ~isempty(findstr([monkey '_' brain_reg '_Units_' MO DA YE '.mat'],horzcat(variables.mat{:})))
    session = cell(12,31);

    WAVES = load([monkey '_' brain_reg '_Units_' MO DA YE]);

    cd(orig_place);

    Wcell = struct2cell(WAVES);
    
    bdfunitinds = vertcat(bdf_in_use.units.id);

    bdf_locs = find(bdfunitinds(:,2) ~= 0 & bdfunitinds(:,2) ~= 255);

    for i = 1:length(Wcell)

        mean_wave = Wcell{i}(4:51);
        std_wave = Wcell{i}(52:end);

        % Check the channel/unit
        chan_unit = Wcell{i}(1:2);  
        if chan_unit == bdf_in_use.units(bdf_locs(i,:)).id    
            bdf_in_use.units(bdf_locs(i,:)).wave = [mean_wave; std_wave];  
        else
            warning('SOMETHING IS WRONG!!');
        end
    end

else
    fprintf('File Not Found\n');
    cd(orig_place);
end



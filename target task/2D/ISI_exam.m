function[perc, acceptable] = ISI_exam(units,ms_cutoff,perc_cutoff,varargin)

if nargin > 1
    perc = zeros(length(units),1);
    for i = 1:length(units)

        spike_times = units{i}(2:end);
        isis = diff(spike_times);

        perc(i) = 100*sum(isis < (ms_cutoff./1000))./length(spike_times);

    end
    acceptable = find(perc < perc_cutoff); 

else
    
    ts = 1:0.1:10;
    perc_profile = zeros(length(units),length(ts));
    for i = 1:length(units)

        spike_times = units{i}(2:end);
        isis = diff(spike_times);

        for tcut = 1:length(ts)
            s_cut = ts(tcut)./1000;
            perc_profile(i,tcut) = 100*sum(isis < s_cut)./length(spike_times);
        end
    end

    isi_under1 = perc_profile < 1;
    
    perc_good_units = 100*sum(isi_under1)./length(units);
    
    figure; plot(ts,perc_good_units);
    xlabel('Temporal cutoff (ms)');
    ylabel('Percentage of Units Acceptable');
    perc = [];
    acceptable = [];
    
end

% %
% [bincounts] = histc(diff(day_units{5})*1000,[0:.25:10]);
% figure; bar([0:0.25:10],bincounts);
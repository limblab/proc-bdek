function[waves] = plot_wave_shapes(bdf,monkey,b_area,mm_dd_yyyy)

orig_place = cd;

if exist(['C:\Users\limblab\Desktop\figures\' monkey '_' b_area '_' mm_dd_yyyy '\Waveforms\'],'file') ~= 7
    mkdir(['C:\Users\limblab\Desktop\figures\' monkey '_' b_area '_' mm_dd_yyyy '\Waveforms\'])
end
cd(['C:\Users\limblab\Desktop\figures\' monkey '_' b_area '_' mm_dd_yyyy '\Waveforms\']);

ids = vertcat(bdf.units.id);

sorted_inds = find(ids(:,2) ~= 0 & ids(:,2) ~= 255);

isitime = 0:0.005:0.150;
plot_isix = (isitime(1:end-1)+isitime(2:end))./2;

waves = zeros(length(sorted_inds),48);
for i = 1:length(sorted_inds)
    
    chan_unit = ids(sorted_inds(i),:);
    
    wave = bdf.units(sorted_inds(i)).wave(1,:);
    waves(i,:) = wave;
    stdwave = bdf.units(sorted_inds(i)).wave(2,:);
    
    mupper = wave + stdwave;
    mlower = wave - stdwave; 
    
    ISI = histc(diff(bdf.units(sorted_inds(i)).ts),isitime);
    ISI(end) = [];
    
    h = figure; subplot(1,2,1); hold on;
    plot(wave,'b'); plot(mupper,'b--'); plot(mlower,'b--');
    ylabel('Millivolts','FontSize',14);
    title(sprintf('Channel %d Unit %d',chan_unit(1),chan_unit(2)),'FontSize',16);
    subplot(1,2,2); bar(plot_isix,ISI); xlabel('ISI (ms)','FontSize',14); ylabel('Count','FontSize',14);
    xlim([0 max(isitime)]);
    
    saveas(h,sprintf('%s_%s_%s_waveform (%d_%d)',monkey,b_area,mm_dd_yyyy,chan_unit(1),chan_unit(2)),'png');
    close(h);
end

cd(orig_place);



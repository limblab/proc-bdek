function track_neurons(session)

alldays = vertcat(session{:});
num_days = length(alldays);

% Loop through days
for i = 1:num_days
    
    ids = vertcat(alldays(i).bdfP.units.id);
    sorted = find(ids(:,2)~= 0 & ids(:,2) ~= 255);
    
    % Loop through sorted units
    for j = 1:length(sorted)
        
        % Obtain and normalize wave shape
        waveform = alldays(i).bdfP.units(sorted(j)).wave(1,:);
        wavestd = alldays(i).bdfP.units(sorted(j)).wave(2,:);
        wavegain = 1/max(waveform); 
        isis = diff(alldays(i).bdfP.units(sorted(j)).ts);
        
        wave.id = ids(sorted(j),:);
        wave.mu = waveform * wavegain;
        wave.std = wavestd * wavegain;
        wave.max = max(waveform);
        wave.min = min(waveform);
        wave.isi = isis(isis < 0.150);
        
        for q = find(1:num_days ~= i)
     
            idsc = vertcat(alldays(q).bdfP.units.id);
            sortedc = find(idsc(:,2)~= 0 & idsc(:,2) ~= 255);
            
            matching_chans = find(ids(sortedc,1)==wave.id(1))
            
            for z = 1:length(matching_chans)
                
                waveformc = alldays(q).bdfP.units(sortedc(matching_chans(z))).wave(1,:);
                wavestdc = alldays(q).bdfP.units(sortedc(matching_chans(z))).wave(2,:);
                wavegainc = 1/max(waveformc); 
                isisc = diff(alldays(q).bdfP.units(sortedc(matching_chans(z))).ts);
                isisc = isisc(isisc < 0.150);
                
                wave_comp.mu = waveformc * wavegainc;
                
                wave_diff = wave_comp.mu - wave.mu;
                
                if sum(wave_diff < abs(wave.std) == length(wave.std))
                % Check ISI
                    
                equiv_ind{i}(j,q) = sortedc(matching_chans(z))
    
            
            
           
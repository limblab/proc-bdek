%% Neuron
neu = 1;
neu_ts = PMd_units{neu};
neu_ts = neu_ts(2:end);

%% Spikes
spikes = zeros(max(round(1000*neu_ts)),1);
spikes(round(1000*neu_ts)) = 1;

trialinds = round(1000*comp_tt(:,[5 7]));

trialset = [];
for i = 1:length(trialinds)
    
    trialset = [trialset trialinds(i,1):trialinds(i,2)];
    
end

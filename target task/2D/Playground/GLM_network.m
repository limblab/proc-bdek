ensemble = alldays(1).PMd_units;
K = 100000;
T = length(alldays(1).kin.pos(:,1));
N = length(ensemble);

B = floor(linspace(1,T,round(T/K)));
for b = 1:length(B)-1
    
    Array = zeros(N,B(b+1)-B(b));
    
    for n = 1:N
        neur = ensemble{n};
        clip = ensemble{n}-1000
        Array(n,
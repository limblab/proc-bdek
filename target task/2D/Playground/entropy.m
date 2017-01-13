ensemble = alldays(1).M1_units;
K = 100000;
T = length(alldays(1).kin.pos(:,1));
N = length(ensemble);

B = 1:100:T;%floor(linspace(1,T,round(T/K)));
for b = 1%:length(B)-1
    
    Array = zeros(N,B(b+1)-B(b));
    TS_TF = alldays(1).kin.pos([B(b),B(b+1)],1)';
    
    for n = 1:N
        neur = ensemble{n}(2:end)-TS_TF(1);
        clip = ceil(1000*neur(neur < diff(TS_TF) & neur > 0));
        Array(n,clip) = 1;
    end
    
    state_i = cell(1,size(Array,2));
    for tb = 1:size(Array,2)
        state_i{tb} = find(Array(:,tb)==1);
    end
    
    % Run comparisons between all 
   
        
        
        
end
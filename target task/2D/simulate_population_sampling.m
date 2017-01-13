% This is a simulation exercise to see if a population sampling theory
% creates population wide neural effects
%% Define some parameters

% Population size
N = 100;    

% Preferred directions
theta0 = 2*pi*rand(N,1);

% Baseline firing rates
base = 15 + 5.*randn(N,1);

% Gains
gain = (0.4 + 0.5*rand(N,1)).*base;
%gain = gamrnd(1.5,1,N,1);

% Number of trials
T1 = 500;   % Condition 1 (center out)
T2 = 500;   % Condition 2 (narrow likelihood)
T3 = 500;   % Condition 3 (broad likelihood)
T = 200;    % Trial duration in milliseconds

% Variance of reach directions
kap1 = 25;
kap2 = 25;

% Simulate reach directions
reach_angle(1:T1) = (pi/4)*floor(8*rand(T1,1))+rand(T1,1)*(10/180*pi);
%reach_angle = 2*pi*rand(T1,1);

for i = T1+(1:T2)
    reach_angle(i) = pi/2 + vonmisrand(kap1);
end
for i = T1+T2+(1:T3)
    reach_angle(i) = pi/2 + vonmisrand(kap2);
end

% How does uncertainty change tuning?
baseshift2 = 0; %(1/5)*gamrnd(5,1,N,1);        % the baseline gets shifted by this factor in Condition 2 (LOW)
baseshift3 = 0; %gamrnd(5,1,N,1);        %0*rand(N,1);      % the baseline gets shifted by this factor in Condition 3 (HIGH)

gainscale2 = 1; %0.1*gamrnd(5,1,N,1);%0.5*rand(N,1);     % the gain gets scaled by this factor in Condition 2 (LOW)
gainscale3 = 1; %0.02*gamrnd(5,1,N,1);%0.1*rand(N,1);         % the gain gets scaled by this factor in Condition 3 (HIGH)

%%
tc1 = gampdf(linspace(0,10,T),2,1) + 0.1; tcourse1 = tc1./max(tc1);
tc2 = exp(linspace(0,1,T)); tcourse2 = tc2./max(tc2);

tcs{1} = tcourse1;
tcs{2} = tcourse2;

tctype = (rand(N,1)<=0.5) + 1;
tgain = [zeros(1,100) linspace(0,1,T-100)];

%% Create center-out trial table
ad(1).tt(:,10) = reach_angle(1:T1);
ad(1).tt(:,3) = 10000;
ad(1).tt(:,5) = 2*(1:T1);
ad(1).tt(:,6) = ad(1).tt(:,5) + T/1000;
ad(1).slices = ad(1).tt(:,10);

ad(2).tt(:,10) = reach_angle((T1+1):end);
ad(2).tt(:,3) = [50*ones(T2,1); 5*ones(T3,1)];
ad(2).tt(:,5) = 2*T1+10 + 2*(1:(T1+T2));
ad(2).tt(:,6) = ad(2).tt(:,5) + T/1000;
ad(2).slices = ad(2).tt(:,10);

%% Method 2
alltts = vertcat(ad(:).tt);
for n=1:N
    %clc; fprintf('Neuron: %03d\n',n);
    
    spikearr = zeros(size(alltts,1),T);
    spiketimes = cell(size(alltts,1),1);
    for tr = 1:size(alltts,1)
        
        clc; fprintf('Neuron: %03d\nTrial: %03d\n',n,tr);
      
        randl = zeros(T,1);
        for zi = 1:T; randl(zi) = vonmisrand(alltts(tr,3)./5); end;
        
        lam = repmat(base(n),T,1) + repmat(gain(n),T,1).*cos((repmat(alltts(tr,10),T,1)+randl) - repmat(theta0(n),T,1));
        
        spikearr(tr, rand(1,T)<=(lam'./1000)) = 1;
        spiketimes{tr} = find(spikearr(tr,:)==1)./1000 + alltts(tr,5);

    end
    
    nspiketimes = horzcat(spiketimes{:});
    ad(1).PMd_units{n} = sortrows(nspiketimes);
    
end


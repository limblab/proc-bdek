% 2014-2-26 
% spike generation with 1 ms refractory period
% reference: M.Okatan, M.A.Wilson, Emery N.Brown, Neural Computation, 2005
% coefficients function:
% a_0 = -2*sin(2*pi*u/0.08)exp(-u/0.04)
% a_+ = 2*sin(2*pi*u/0.06)exp(-u/0.04)
% a_- = -3*sin(2*pi*u/0.12)exp(-u/0.04)

clc;clear all;close all;
M = 120;% history effect   120ms
C = 4;%the number of neurons
a0 = ones(4,1)*log(5);% spontaneous firing rate 5Hz

T = 100*1000;% length of the data 60s
dt = 1;% 1ms

refrac_flag1 = false;
refrac_flag2 = false;
refrac_flag3 = false;
refrac_flag4 = false;

spike_time = zeros(C,T+1);
spike_rate = zeros(C,T+1);
% lambda = zeros(C,T);
spike_time(:,1) = 0.001*rand(C,1);% at time 0
spike_rate(:,1) = 1;% I_0=1;
for ts = 1:dt:T
    
    if  ts <= M      
        u = (0:dt:ts-1)/1000;%in seconds
        W = (0:dt:ts-1)+1;
    else
        u = (1:dt:M)/1000;
        W = (ts-M+1:dt:ts);
    end
    a_0 = -2*sin(2*pi*u'/0.08).*exp(-u'/0.04);%itself
    a_e = 2*sin(2*pi*u'/0.06).*exp(-u'/0.04);%excitation
    a_i = -3*sin(2*pi*u'/0.12).*exp(-u'/0.04);%inhibition
    a_w = zeros(size(a_0));% the neuron without connections
    a1 = [a_0;a_i;a_w;a_e];
    a2 = [a_e;a_0;a_i;a_w];
    a3 = [a_w;a_e;a_0;a_i];
    a4 = [a_i;a_w;a_e;a_0];
    
    I = [spike_rate(1,W),spike_rate(2,W),spike_rate(3,W),spike_rate(4,W)];
    lambda1 = exp(a0+I*a1);
    lambda2 = exp(a0+I*a2);
    lambda3 = exp(a0+I*a3);
    lambda4 = exp(a0+I*a4);
     
    t = ts + 1;
    %%%%%%%%%neuron 1%%%%%%%%%
    if refrac_flag1
        lambda1 = 0;
        refrac_flag1 = false;
    end    
    if lambda1*dt/1000 > rand
        spike_time(1,t) = ts;
        spike_rate(1,t) = 1;
        refrac_flag1 = true;
    end
    
    
    %%%%%%%%%neuron 2%%%%%%%%%
    if refrac_flag2
        lambda2 = 0;
        refrac_flag2 = false;
    end    
    if lambda2*dt/1000 > rand
        spike_time(2,t) = ts;
        spike_rate(2,t) = 1;
        refrac_flag2 = true;
    end
    
    
    %%%%%%%%%neuron 3%%%%%%%%%
    if refrac_flag3
        lambda3 = 0;
        refrac_flag3 = false;
    end    
    if lambda3*dt/1000 > rand
        spike_time(3,t) = ts;
        spike_rate(3,t) = 1;
        refrac_flag3 = true;
    end
    
    
    %%%%%%%%%neuron 4%%%%%%%%%
    if refrac_flag4
        lambda4 = 0;
        refrac_flag4 = false;
    end
    if lambda4*dt/1000 > rand
        spike_time(4,t) = ts;
        spike_rate(4,t) = 1;
        refrac_flag4 = true;
    end
    
end

% for i = 1:C
%     
%     Idx = find(spike_time(i,:));
%     isi = diff(spike_time(i,Idx));
%     maxIsi = ceil(max(isi) / 1000) * 1000;
%     bins = -0.5 : maxIsi - 0.5;
%     subplot(2,2,i);hist(isi, bins);
%     isiHist = hist(isi, bins);
%     isiHist = isiHist / sum(isiHist);
%     
%     %we need the cumultative sum for the survivor function
%     cumTih = cumsum(isiHist);
%     survivorFunc = 1 - cumTih;
%     
%     %hazard functions is the ISI divided by the survivor function
%     hazard = isiHist ./ survivorFunc;
% 
% end

spk_times = cell(1,1);
for i = 1:C
    idx = find(spike_time(i,:));
    spk_times{i} = spike_time(i,idx)/1000; % in seconds
end

% the number of spikes
length(find(spike_rate(1,:)))
length(find(spike_rate(2,:)))
length(find(spike_rate(3,:)))
length(find(spike_rate(4,:)))

save('syntheticSpk','spk_times');

% 2014-2-26 
% glmfit to analyze the functional connectivity of the four synthetic
% neurons
clc;clear all;close all;

load syntheticSpk;

Tmax = 100;
delta = 0.01;% binsize = 10ms
times = 0:delta:Tmax;
times = times';
bin_times = [-.5*delta; times+.5*delta];
N = length(bin_times)-1;
C = size(spk_times,2);
spk_rates = zeros(N,C);

for i = 1:C
    spk_times_c = spk_times{i};
    for j = 1:N
        if j==N
            spk_rates(j,i) = sum(spk_times_c>=bin_times(j)&spk_times_c<=bin_times(j+1));
        else
            spk_rates(j,i) = sum(spk_times_c>=bin_times(j)&spk_times_c<bin_times(j+1));
        end
    end
end

save('spikeRate001.mat','spk_rates');

M = 12;
duration = 51:size(spk_rates,1)-M;

train_len = length(duration);
X = zeros(train_len,M*C);
for l = 1:train_len
    tmp = [];
    for n = 1:C
        tmp = [tmp;spk_rates(duration(l):duration(l)+(M-1),n)];
    end
    X(l,:) = tmp';
end

distr = 'poisson';
link = 'log';
output = cell(1,1);
for n = 1:C
   y = spk_rates(duration+M,n);
   [b, dev, stats] = glmfit(X,y,distr,'link',link,'constant','off');
   b = real(b);    
   output{n}.b = reshape(b,M,C);
end

for n = 1:C
    subplot(C,C,(n-1)*4+1);
    plot(output{n}.b(:,1));
    xlim([0 12]);
    subplot(C,C,(n-1)*4+2);
    plot(output{n}.b(:,2));
    xlim([0 12]);
    subplot(C,C,(n-1)*4+3);
    plot(output{n}.b(:,3));
    xlim([0 12]);
    subplot(C,C,(n-1)*4+4);
    plot(output{n}.b(:,4));
    xlim([0 12]);
end

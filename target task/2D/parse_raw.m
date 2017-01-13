loc = 'C:\Users\limblab\Desktop\Raw Data\';
base_name = 'Mihili_PMD_05202014_UNT2D_001';

spike = cell(96,1);
for ch = 1:96
    clc;fprintf('channel: %d/96\n',ch);
    clearvars -except spike loc base_name ch
    dat = importdata([loc base_name sprintf('_%03d.txt',ch)]);
   
    spike{ch,1} = dat.data(:,1);
    spike{ch,2} = dat.data(:,2:end);
end

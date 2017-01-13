function[cont_signal,filtmax] = train2cont(train,maxISI,causal,varargin)

% train2cont(train,reach) converts a spike train (raster) to a continuous
% signal using Gaussian convolution. The rise time will be approximately
% 1.5*maxISI

% train: Input spike train (zeros and ones)
% maxISI: represents two standard deviations of the Gaussian kernel (in ms) 
%        Example: maxISI = 50 means that the amplitude of the continuous
%        signal will only significantly increase when an ISI is less than
%        50 ms. 
% causal: If (1), then maxISI is doubled and only the causal half of the
% gaussian is used. 

sd = maxISI/2;
N = 6*sd+1;
if nargin > 2 && causal
    
    causfilt = fspecial('gaussian',[2*N 1],2*sd)*2;
    causfilt(1:floor(end/2)) = 0; 

    cont_signal = conv(train,causfilt,'same')*1000;
    filtmax = max(causfilt)*1000;
else
    gaussfilt = fspecial('gaussian',[2*N 1],2*sd);
    cont_signal = conv(train,gaussfilt,'same')*1000;
    filtmax = max(gaussfilt)*1000;
end


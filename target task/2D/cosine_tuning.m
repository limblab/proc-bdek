function[B,fitp,pref_d,offset,magnitude] = cosine_tuning(tt,train,t1,t2,noise_dist,alignment)

if strcmp(alignment,'target');
    align_col = 5;
elseif strcmp(alignment,'go');
    align_col = 6;
else
    align_col = 5;
    fprintf('alignment not recognized: target used\n');
end

% Find directions for use in tuning 
direcs = tt(:,10);

% Calculate binned
[rast,y] = bin_rast(tt,train,t1,t2,1,align_col);
y = y./((t2-t1)/1000);

X = [cos(direcs) sin(direcs)];

[B,dev,stats] = glmfit(X,y,noise_dist);

offset = B(1);
pref_d = atan(B(3)/B(2));
magnitude = B(2)/cos(pref_d);

stats.p;

% Test against constant model
[Bc,devc] = glmfit(ones(length(direcs),1),y,noise_dist,'constant','off');
fitp = 1 - chi2cdf(devc - dev,2);

if strcmp(noise_dist,'poisson')
    noise_distval = 'log';
elseif strcmp(noise_dist,'normal')
    noise_distval = 'identity';
else
    fprintf('Type not found\n');
    noise_distval = 'identity';
end

% xs = 0:0.01:2*pi;
% figure; plot(direcs,y,'b.');
% hold on;
% plot(xs,glmval(B,[cos(xs)' sin(xs)'],noise_distval),'r');

end





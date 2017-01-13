if 1
l1 = 20; %35;
l2 = 25; %30;

theta1 = 0:0.01:pi/2;
theta2 = 0:0.01:pi;

zeroangs = [1.28 1.57];

[THETA1,THETA2] = meshgrid(theta1,theta2);
X = l1*cos(THETA1) + l2*cos(THETA1+THETA2);
Y = l1*sin(THETA1) + l2*sin(THETA1+THETA2);

data1 = [X(:) Y(:) THETA1(:)];
data2 = [X(:) Y(:) THETA2(:)];

zeropoint = data1(data1(:,3)==zeroangs(1) & data2(:,3)==zeroangs(2),1:2);

%%
angs = .001:0.001:2*pi;
[TS,Delta_ts] = deal(zeros(length(angs),2));
for i = 1:length(angs)
    
    XY = zeropoint + [cos(angs(i)) sin(angs(i))];
    
    X = XY(1); Y = XY(2);
    c2 = (X.^2 + Y.^2 - l1^2 - l2^2)/(2*l1*l2);
    s2 = sqrt(1-c2.^2);
    THETA2D = atan2(s2,c2);

    k1 = l1 + l2.*c2;
    k2 = l2*s2;
    THETA1D = atan2(Y,X) - atan2(k2,k1);
    
    TS(i,:) = [THETA1D, THETA2D]; 
    Delta_ts(i,:) = TS(i,:)-zeroangs;
end
dir_ang = atan2(Delta_ts(:,2),Delta_ts(:,1));
[yp,xp] = hist(dir_ang,100);
prior_f = @(x) interp1(xp,yp,x,'spline');

%%
K_l = 5; 
ts = -pi:0.001:pi;

% Prior distribution
p_dist = prior_f(ts);

[p_of_chosen,model_map] = deal(zeros(length(Bdir),1));
for i = 1:length(Bdir) % Loop through trials
    clc; fprintf('Trial: %d/%d\n',i,length(Bdir));
    [theta_hat,MAP] = deal(zeros(1000,1));
    for smp = 1:length(theta_hat) % Loop through monte carlo samples
        
        theta_hat(smp) = Bdir(i) + vonmisrand(K_l); % noisy estimated samples of theta
        l_dist = circ_vmpdf(ts,theta_hat(smp),K_l); % Likelihood function for given theta_hat
        
        posterior = l_dist.*p_dist'; % Posterior distribution
        
        MAP(smp) = ts(posterior==max(posterior));
        
    end
    
    p_of_chosen(i) = circ_vmpdf(Rdir(i),circ_mean(MAP),circ_kappa(MAP));
    model_map(i) = circ_mean(MAP);
 
end
   
end
%%
tic; 
K_l = 5; 
d_theta = 0.01;
reps = 200;

dtc = 1/d_theta;
ts = -pi:d_theta:pi;
% Prior distribution
p_dist = prior_f(ts);
pdmat = repmat(p_dist,reps,1);

% Likelihood distribution
l_dist = circ_vmpdf(ts,0,K_l);
ldmat = zeros(length(ts));
for i = 1:length(ts)
    ldmat(i,:) = circshift(l_dist,i-round(length(ts)/2))';
end
% Preconstruct von mises lookup table
vm_lookup = zeros(1000,1);
for i = 1:1000
    vm_lookup(i) = round(dtc*vonmisrand(K_l))/dtc;
end

[p_of_chosen,model_map] = deal(zeros(length(Bdir),1));
theta2ind = @(x) mod(floor(dtc*x + length(ts)/2),length(ts)-1)+1;
indexing_mat = repmat(1:length(ts),reps,1);
for i = 1:length(Bdir) % Loop through trials
    %clc; fprintf('Trial: %d/%d\n',i,length(Bdir));

    theta_hats = repmat(floor(dtc*Bdir(i))/dtc,reps,1) + datasample(vm_lookup,reps); % noisy estimated samples of theta
    theta_hat_inds = theta2ind(theta_hats);

    posteriors = ldmat(theta_hat_inds,:).*pdmat;    
    MAP = ts(indexing_mat(posteriors == repmat(max(posteriors,[],2),1,size(posteriors,2))))';
      
    %p_of_chosen(i) = circ_vmpdf(Rdir(i),circ_mean(MAP),circ_kappa(MAP));
    model_map(i) = circ_mean(MAP);

end

% theta_hats_inds = theta2ind(repmat((floor(1000*Bdir)/1000)',1000,1) + reshape(datasample(vm_lookup,1000*length(Bdir)),1000,[]));
% theta_hat_vec = reshape(theta_hats_inds,[],1);
% like_vec = ldmat(theta_hat_vec,:);
% like_mat = reshape(like_vec,1000,[]);

te = toc; fprintf('Time: %d m %d s\n',floor(te/60),round(te-60*floor(te/60)));


%% Fit likelihood kappa to behavioral data



%fit_kap = fminsearch(helper_func,1);

[~,samps] = bootstrp(20,@(x) x,1:length(Bdir));
fits_kap = zeros(size(samps,2),1);
[es,mods] = deal(cell(size(samps,2),1));
ers = zeros(size(samps,2),20);
for i = 1:size(samps,2)
    clc; fprintf('%d/%d\n',i,size(samps,2));
    Bdir_samp = Bdir(samps(:,i));
    Rdir_samp = Rdir(samps(:,i));
    
%     helper_func = @(K) Bayes_fitting_func(Bdir_samp,Rdir_samp,prior_f,K);
%     %fits_kap(i) = fminunc(helper_func,1,optimset('display','iter','MaxIter',10,'DiffMinChange',0.1,'DiffMaxChange',1));
%     fits_kap(i) = fmincon(helper_func,1,[],[],[],[],0.1,20,[],optimset('algorithm','interior-point',...
%         'display','off','MaxIter',10,'DiffMinChange',0.1,'DiffMaxChange',1));
%     [es{i},mods{i}] = Bayes_fitting_func(Bdir_samp,Rdir_samp,prior_f,fits_kap(i));
    
    test_kaps = 1:1:20;
    for j = 1:length(test_kaps)
        clc; fprintf('%d/%d - %d/%d\n',i,size(samps,2),j,20);
        [ers(i,j),mods{i}(:,j)] = Bayes_fitting_func(Bdir_samp,Rdir_samp,prior_f,test_kaps(j));
    end
end


%fit_kap2 = fminunc(helper_func,1,optimset('display','iter','MaxIter',10,'DiffMinChange',0.1,'DiffMaxChange',1));
%%
[~,model_map] = Bayes_fitting_func(Bdir,Rdir,prior_f,7);

figure; hold on; plot(Bdir,model_map,'.'); plot(Bdir,Rdir,'r.');

%%
test_kaps = 1:1:20;
for i = 1:length(test_kaps)
    clc; fprintf('%d/%d\n',i,length(test_kaps));
    [sse(i),mods{i}] = Bayes_fitting_func(Bdir,Rdir,prior_f,test_kaps(i));
end

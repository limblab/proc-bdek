kp = 3; % Prior std

% Feedback kappas
like_kappas = [25 2];

ntrials = 1000; %Trials per session     
sigp_hats = 0:0.1:4; %Range of prior estimates (strategies)
strategies = [3]; %

%% Mr. F (the fake monkey)
centroid_estimation_error = 0; % percent var
motor_error = 0; % std (constant)


for i = 1:length(ntrials)
    
    tloc(i) = vonmisrand(kp);

    for j = 1:length(like_kappas)
        feedloc{j}(i) = tloc(i) + vonmisrand(like_kappas(j));
        
        opt_weight = 1/(kp/like_kappas(j) + 1); 

    
    end
end


%%
muls = (pi/2):0.01:pi;
kls = 0.5:0.01:2;
[d,d2] = deal(zeros(length(muls),length(kls)));
for i = 1:length(muls)

    mu_l(i) = muls(i);
    
    for j = 1:length(kls)
       
        d(i,j) = dfunc(mu_p,mu_l(i),k_p,kls(j)*25);
        d2(i,j) = dfunc(mu_p,mu_l(i),k_p,kls(j)*250);
    end
end
    
tarray = repmat(muls',1,length(kls));

errs1 = sum(d<(5/180*pi),2);
errs2 = sum(d2<(5/180*pi),2);
    
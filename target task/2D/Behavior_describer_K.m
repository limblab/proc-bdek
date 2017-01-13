tt = alldays(2).tt;

[delta_k, best_k, new_B_est] = deal(zeros(size(tt,1),1));
for i = 1:size(tt,1)

    movedir = tt(i,10);    
    %%%%%%%%%%%%%%%%%%% LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lines = alldays(2).slices(i,:);
    n = length(lines);
    
    % CALCULATE PARAMETERS OF LIKELIHOOD VON MISES
    % From "Statistical Analysis of Circular Data" N.I. Fisher (p 88-9)
    R = norm([sum(cos(lines)) sum(sin(lines))]); % Resultant vector
    Rbar = R/n; %  mean Resultant vector (4.38)
    
    % Calculate k_ml according to Rbar (4.40)
    if Rbar < 0.53
        k_ml = 2*Rbar + Rbar.^3 + 5*Rbar^5/6;
    elseif (Rbar >= 0.53 && Rbar <= 0.85)
        k_ml = -0.4 + 1.39*Rbar + 0.43/(1-Rbar);
    elseif Rbar > 0.85
        k_ml = 1/(Rbar^3 - 4*Rbar^2 + 3*Rbar);
    else
        fprintf('Problem with k_ml (likelihood) calculation: Check cases\n');
    end
    
    % Adjust khat if n <= 15 (4.41)
    if n <= 15
        if k_ml < 2
            khat = max([k_ml - 2*(n*k_ml)^(-1)  ,  0]);
        elseif k_ml >= 2
            khat = (n-1)^3*k_ml/(n^3 + n); 
        else
            fprintf('Problem with khat (likelihood) calculation: Check cases\n');
        end
    else
        khat = k_ml;
    end
    
    k_likelihood = sqrt(Rbar*khat); % (4.42)
    mu_likelihood = circ_mean(lines');
    
    %%%%%%%%%%%%%%%%%%%%%% PRIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prev_locs = tt(:,2);
    n = length(prev_locs);

    % CALCULATE PARAMETERS OF PRIOR VON MISES
    % From "Statistical Analysis of Circular Data" N.I. Fisher (p 88-9)
    R = norm([sum(cos(prev_locs)) sum(sin(prev_locs))]); % Resultant vector
    Rbar = R/n; %  mean Resultant vector (4.38)

    % Calculate k_ml according to Rbar (4.40)
    if Rbar < 0.53
        k_ml = 2*Rbar + Rbar.^3 + 5*Rbar^5/6;
    elseif (Rbar >= 0.53 && Rbar <= 0.85)
        k_ml = -0.4 + 1.39*Rbar + 0.43/(1-Rbar);
    elseif Rbar > 0.85
        k_ml = 1/(Rbar^3 - 4*Rbar^2 + 3*Rbar);
    elseif isnan(Rbar)
        k_ml = NaN;
    else    
        fprintf('Problem with k_ml (prior) calculation: Check cases\n');
    end

    % Adjust khat if n <= 15 (4.41)
    if n <= 15
        if k_ml < 2
            khat = max([k_ml - 2*(n*k_ml)^(-1)  ,  0]);
        elseif k_ml >= 2
            khat = (n-1)^3*k_ml/(n^3 + n); 
        elseif isnan(k_ml)
            khat = NaN;
        else
            fprintf('Problem with khat (prior) calculation: Check cases\n');
        end
    else
        khat = k_ml;
    end

    k_prior = sqrt(Rbar*khat); % (4.42)
    mu_prior = circ_mean(prev_locs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    kfunc = @(l_k) (circ_dist(Bayes_circ_trial(mu_likelihood,l_k,mu_prior,k_prior),movedir)).^2;
    
    best_k(i) = fminbnd(kfunc,-1000,1000);
    delta_k(i) = best_k(i) - k_likelihood;
    
    new_B_est(i) = Bayes_circ_trial(mu_likelihood,best_k(i),mu_prior,k_prior);
    
    clc; fprintf('trial: %d\n',i);

end
    
    
    
function[k] = calc_kappa(target_angs)

lines = target_angs;
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

k = khat;

end
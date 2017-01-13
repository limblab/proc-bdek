function vmvec = vonmisrand_vec(n,K)

% if length(K) == 1:
%       generate vector of length n that follows von Mises distribution 
%       with mean 0 and kappa K
% if length(K) > 1:
%       generate vector of length K where the nth element of vmvec is a von
%       Mises random sample with mean 0 and kappa K(n)

if length(K) == 1
    vmvec = zeros(n,1);
    for i = 1:n
        vmvec(i) = vonmisrand(K);
    end

else
    vmvec = zeros(length(K),1);
    
    for i = 1:length(K)
        vmvec(i) = vonmisrand(K(i));
    end
end

end
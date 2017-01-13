function [max_prob] = prob_help_func(x,targ_means,num_targs)

pAgivenTs = poisspdf(repmat(x,num_targs,1),targ_means);
pAgivenTs(pAgivenTs==0) = eps;
marg = nansum(pAgivenTs,1);
pTgA = pAgivenTs./repmat(marg,num_targs,1);

lograts = sum(log(pTgA),2);

max_prob = max(lograts);
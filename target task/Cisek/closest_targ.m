function[total_p] = closest_targ(slide,foraverage,againstaverage,x)

valids = find(prod([foraverage againstaverage slide],2)~=0);

Spike = slide(valids);
Fspike = foraverage(valids);
Aspike = againstaverage(valids);

nearprob_func = @(x) max([poisspdf(Spike,x(1)*Fspike+x(2)),...
                          poisspdf(Spike,x(1)*Aspike+x(2))],[],2);
             
total_p = -mean(nearprob_func(x));

 
function [ys1, ys2] = plottick(x,y,numtrials,spread,cond,color)
    
    if numtrials==spread
        plot(x,y,'.','Color',color','MarkerSize',50);
    else
        range = linspace(0,spread,numtrials+1);
        ys1 = spread*(cond-1)+range(mod(y,spread));
        ys2 = spread*(cond-1)+range(mod(y,spread)+1);
        
        plot([x;x],[ys1;ys2],'Color',color,'LineWidth',1);
    end
end
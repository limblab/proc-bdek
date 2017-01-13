function[xp,yp] = linehist(x,y,color2plot,varargin)


if numel(y)>1
    
    n = histc(y,x);
    b = x;
    
    d = mode(diff(b));
    drep = [-(d/2)*ones(1,length(b)); (d/2)*ones(1,length(b))];

    xp = reshape(repmat(reshape(b,1,[]),2,1)+drep,1,[]) ;
    yp = [0 reshape(repmat(reshape(n,1,[]),2,1),1,[]) 0];
    
    xp = [xp(1) xp(1:end) xp(end)];
    
    if nargin > 2
        plot(xp,yp,color2plot,'LineWidth',2);
    else
        plot(xp,yp,'LineWidth',2);
    end
else
    
    [n,b] = hist(x,y);

    d = mode(diff(b));
    drep = [-(d/2)*ones(1,length(b)); (d/2)*ones(1,length(b))];
    
    xp = reshape(repmat(reshape(b,1,[]),2,1)+drep,1,[]) ;
    yp = [0 reshape(repmat(reshape(n,1,[]),2,1),1,[]) 0];
    
    xp = [xp(1) xp(1:end) xp(end)];
    
    if nargin > 2
        plot(xp,yp,color2plot,'LineWidth',2);
    else
        plot(xp,yp,'LineWidth',2);
    end
end
    
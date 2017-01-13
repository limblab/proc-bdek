function [speeds,xp,yp,xvel,yvel,xa,ya] = kin_exam(bdf,tt)

%movetimes = tt(:,6)-tt(:,5);
movetimes = 1000;% tt(:,6) + 1000;
%maxmove = ceil(1000*max(movetimes));
maxmove = 1000;

[xvel,yvel,xp,yp,xa,ya,speeds] = deal(nan(length(tt),maxmove));

for i = 1:length(tt)
    
    timego = tt(i,5);
    %timeend = tt(i,6);
    %timeend = tt(i,6)+1;
    
    bdfgo = find(bdf.pos(:,1)>timego,1,'first');
    bdfend = bdfgo+1400;
    %bdfend = find(bdf.pos(:,1)>timeend,1,'first');
    
    xvel(i,1:bdfend-bdfgo+1) = bdf.vel(bdfgo:bdfend,2)';
    yvel(i,1:bdfend-bdfgo+1) = bdf.vel(bdfgo:bdfend,3)';
    xp(i,1:bdfend-bdfgo+1) = bdf.pos(bdfgo:bdfend,2)';
    yp(i,1:bdfend-bdfgo+1) = bdf.pos(bdfgo:bdfend,3)';
    xa(i,1:bdfend-bdfgo+1) = bdf.acc(bdfgo:bdfend,2)';
    ya(i,1:bdfend-bdfgo+1) = bdf.acc(bdfgo:bdfend,3)';
    
    speeds(i,1:bdfend-bdfgo+1) = sqrt(bdf.vel(bdfgo:bdfend,2).^2 + bdf.vel(bdfgo:bdfend,3).^2)';
%     xpos(i,1:bdfend-bdfgo+1) = BDF.pos(bdfgo:bdfend,2)';
%     ypos(i,1:bdfend-bdfgo+1) = BDF.pos(bdfgo:bdfend,3)';
    
end
    
    
    

    
function [speeds,xp,yp,xvel,yvel,xa,ya] = kin_exam(bdf,tt,move_1,move_2,align)

if length(align) > 1

    [xvel,yvel,xp,yp,xa,ya,speeds] = deal(cell(size(tt,1),1));
    for i = 1:size(tt,1)

        timego = tt(i,align(1));
        timeend = tt(i,align(2));

        bdfgo = find(bdf.pos(:,1)>(timego),1,'first');
        bdfend = find(bdf.pos(:,1)<(timeend),1,'last');

        xvel{i} = bdf.vel(bdfgo:bdfend,2)';
        yvel{i} = bdf.vel(bdfgo:bdfend,3)';
        xp{i} = bdf.pos(bdfgo:bdfend,2)';
        yp{i} = bdf.pos(bdfgo:bdfend,3)';
        xa{i} = bdf.acc(bdfgo:bdfend,2)';
        ya{i} = bdf.acc(bdfgo:bdfend,3)';

        speeds{i} = sqrt(bdf.vel(bdfgo:bdfend,2).^2 + bdf.vel(bdfgo:bdfend,3).^2)';

    end
else
    
    %movetimes = tt(:,6)-tt(:,5);
    %movetimes = 600;% tt(:,6) + 1000;
    %maxmove = ceil(1000*max(movetimes));
    maxmove = move_2 - move_1;
    [xvel,yvel,xp,yp,xa,ya,speeds] = deal(nan(length(tt),maxmove));
    for i = 1:length(tt)

        timego = tt(i,align);
        %timeend = tt(i,6);
        %timeend = tt(i,6)+1;

        bdfgo = find(bdf.pos(:,1)>(timego+move_1./1000),1,'first');
        bdfend = bdfgo+maxmove;
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
end
    
    
    

    
function[p] = ang2poly(angle,width,r1,r2)

    xs = [r1.*cos((angle-width/2):.01:(angle+width/2)) r2.*cos((angle+width/2):-.01:(angle-width/2))]';
    ys = [r1.*sin((angle-width/2):.01:(angle+width/2)) r2.*sin((angle+width/2):-.01:(angle-width/2))]';
    
    p = [xs ys zeros(length(xs),1)];
end

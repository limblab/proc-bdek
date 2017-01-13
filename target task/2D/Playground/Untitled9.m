x = 3;

r = 7;

angls = 0:360;
[mE,mL,mLR,mA,mER] = deal(zeros(length(angls),1));
for i = 1:length(angls)
    a = deg2rad(180-angls(i));
    mL(i) = sqrt(r^2 + x^2 - 2*r*x*cos(a));
    mLR(i) = sqrt(r^2 + x^2 - 2*r*x*cos(pi-a));
    
    mE(i) = rad2deg(asin(x*sin(a)/mL(i)));
    mER(i) = -rad2deg(asin(x*sin(pi-a)/mLR(i)));
    
    mA(i) = rad2deg(pi+angls(i)-mE(i));%asin(r*sin(a)/mL(i));

end

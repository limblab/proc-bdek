
l1 = 20; %35;
l2 = 25; %30;

theta1 = 0:0.01:pi/2;
theta2 = 0:0.01:pi;

zeroangs = [1.06 0.77]; %[1.28 1.57];

[THETA1,THETA2] = meshgrid(theta1,theta2);

X = l1*cos(THETA1) + l2*cos(THETA1+THETA2);
Y = l1*sin(THETA1) + l2*sin(THETA1+THETA2);

data1 = [X(:) Y(:) THETA1(:)];
data2 = [X(:) Y(:) THETA2(:)];

zeropoint = data1(data1(:,3)==zeroangs(1) & data2(:,3)==zeroangs(2),1:2);


angs = .01:0.01:2*pi;
[TS,posits] = deal(zeros(length(angs),2));
for i = 1:length(angs)
    
    XY = zeropoint + [cos(angs(i)) sin(angs(i))];
    
    posits(i,:) = XY;
    
    X = XY(1);
    Y = XY(2);

    c2 = (X.^2 + Y.^2 - l1^2 - l2^2)/(2*l1*l2);
    s2 = sqrt(1-c2.^2);
    THETA2D = atan2(s2,c2);

    k1 = l1 + l2.*c2;
    k2 = l2*s2;
    THETA1D = atan2(Y,X) - atan2(k2,k1);
    
    TS(i,:) = [THETA1D, THETA2D];
    
end
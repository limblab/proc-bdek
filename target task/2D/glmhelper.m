function dir = glmhelper(x)
b = glmfit([cos(x(:,1)),sin(x(:,1))],x(:,2));
dir = atan2(b(3),b(2));
end
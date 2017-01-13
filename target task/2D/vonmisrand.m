function[s] = vonmisrand(K)

tau = 1+sqrt(1+4*K^2);
rho = (tau-sqrt(2*tau))/(2*K);
r   = (1+rho^2)/(2*rho);

accept = 0;
num_iter=0;
fail_iter = 50;
while ((accept==0)&&(num_iter < fail_iter)) 
		u2 = rand;
        
        z = cos(rand*pi);
		f = (1+r*z)/(r+z);
		c = K*(r-f);
       
        
		if (c*(2-c)-u2)>0
			accept = 1;
        elseif ((log(c/u2)+1-c)>=0) 
			accept = 1;
        else
			num_iter = num_iter+1;
        end
end
		
u3 = rand-0.5;
if (u3 < eps) 
    s = (-acos(f)); 
else 
    s = (sign(u3))*acos(f); 
end

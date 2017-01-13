function[d_ang] = angle_diff(ang1,ang2)
% d_ang = angle_diff(ang1,ang2) gives the angle FROM ang1 TO ang2.
% 
% If computing errors between an angle A and reference angle R, use
% d_ang = angle_diff(R,A); This way, negative distances represent an
% UNDERSHOOT of the reference angle. 

if size(ang1,2)>size(ang1,1)
    ang1 = ang1';
end
if size(ang2,2)>size(ang2,1)
    ang2 = ang2';
end

p1 = [cos(ang1) sin(ang1)];
p2 = [cos(ang2) sin(ang2)];

dotp = p1(:,1).*p2(:,1)+p1(:,2).*p2(:,2); % Dot product
crossp = p1(:,1).*p2(:,2)-p2(:,1).*p1(:,2); % Cross product

d_ang = acos(dotp).*sign(crossp); % Find angle and apply correct sign

end


function[centroid] = mean_angle(angles,angle_type)
% centroid = mean_angle(angles,angle_type) calculates the
% centroid angle of the set of input angles. 
% 
% angle_type can be either 'rads' or 'degs' and will return the centroid 
% in the same units

%% Check angle type (radians or degrees) and convert/or return error
if strcmp(angle_type,'degs');
    angles = angles.*pi./180;
elseif ~strcmp(angle_type,'rads') % Return error if improper angle type
    error('Unknown angle type "%s"',angle_type);
end

%% Calculate centroid
av_ang = atan2(mean(sin(angles)),mean(cos(angles))); %Calculate avg vector

%% Output centroid and convert units if necessary
if strcmp(angle_type,'degs')
    centroid = av_ang(1).*180./pi; %Convert output if 'degrees'
else
    centroid = av_ang(1);
end

end

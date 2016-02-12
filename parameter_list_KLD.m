function[param_list] = parameter_list_KLD()

%parameter_list_KLD creates a struct containing input parameters for
%KLD_nongauss.m



param_list = struct('x_pos', 'omit','y_pos', 'omit', ... 
                    'x_vel', 'include', 'y_vel','include', ...
                    'x_acc', 'omit', 'y_acc', 'omit', ...
                    'x_force', 'omit', 'y_force', 'omit', ...
                    'samples', 0.2,'window', 3000,...
                    'unit', 1,'block_param', 0.10, ...
                    'exclusion_type','axes','resp_param',0.1,...
                    'resp_axes','omit','graph', 'yes','show_components','no');
end


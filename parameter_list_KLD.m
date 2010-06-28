function[param_list] = parameter_list_KLD()

param_list = struct('x_pos', 'include','y_pos', 'include', ... 
                    'x_vel', 'include', 'y_vel','include', ...
                    'x_acc', 'omit', 'y_acc', 'omit', ...
                    'x_force', 'omit', 'y_force', 'omit', ...
                    'samples', 1000,'window', 4000, 'unit', 1, ...
                    'block_param', 0.05, 'graph', 'yes');
end


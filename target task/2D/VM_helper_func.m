function[curve] = VM_helper_func(xsys)

[~,~,curve] = VM_fit3(xsys(:,1),xsys(:,2),[-pi:0.01:pi]);

end
function[fwhm] = VM_helper_func_fwhm(xsys)

[~,~,~,~,~,fwhm] = VM_fit(xsys(:,1),xsys(:,2),[-pi:0.01:pi]);

end
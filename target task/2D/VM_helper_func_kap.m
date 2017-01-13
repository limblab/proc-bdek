function[ps] = VM_helper_func_kap(xsys)

[~,~,~,base,gain,fwhm] = VM_fit3(xsys(:,1),xsys(:,2),-pi:0.01:pi);
ps = [base gain fwhm];

end
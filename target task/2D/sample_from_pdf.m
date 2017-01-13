function[samples] = sample_from_pdf(thetas,pdf,n)

pdfr = pdf + .0001*rand(size(pdf));
pdfn = pdfr./sum(pdfr);
%pdfne = pdf./sum(pdf) + eps*rand(size(pdf));
cdf = sum(tril(repmat(pdfn,length(pdfn),1)),2);

unrand = rand(n,1);

samples = interp1(cdf,thetas,unrand,'linear','extrap');

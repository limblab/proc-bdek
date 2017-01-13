function[cdf] = pdf2cdf(pdf)

pdfn = (pdf - min(pdf))./(max(pdf)-min(pdf));

cdf = sum(tril(repmat(pdfn,length(pdfn),1)),2);
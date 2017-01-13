function[curve, pfits, base, gain, fwhm] = extract_boothelper(outstructL,outstructH)

[curve,pfits,base,gain,fwhm] = deal(cell(length(outstructL),1));
for i = 1:length(outstructL)
    
    curve{i}{1} = outstructL{i}(1:629,:);
    curve{i}{2} = outstructH{i}(1:629,:);
    pfits{i}{1} = outstructL{i}(630:633,:);
    pfits{i}{2} = outstructH{i}(630:633,:);
    
    base{i} = [outstructL{i}(634,:);outstructH{i}(634,:)];
    gain{i} = [outstructL{i}(635,:);outstructH{i}(635,:)];
    fwhm{i} = [outstructL{i}(636,:);outstructH{i}(636,:)];
end
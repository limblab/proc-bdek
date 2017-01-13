C = nchoosek(1:11,6);
rinds = randperm(size(C,1));
reps = size(C,1);

ccf_dRES = zeros(reps,4);
[dRESslopes, dslopesslopes] = deal(zeros(reps,4));
ccf_dslopes = zeros(reps,3);
for bsamp = 1:reps;
    
%     G = randperm(23,15);
    G = C(rinds(bsamp),:);

    %% Do
    for tp = 1:length(AV)


    av_PD = AV{tp}{1};
    av_OD = AV{tp}{2};
    av_ORTH = AV{tp}{3};

    d_PD = EF{tp}{1}(:,3);
    d_OD = EF{tp}{2}(:,3);
    d_ORTH = EF{tp}{3}(:,3);

    bnd_PD = EF{tp}{1}(:,1:2);
    bnd_OD = EF{tp}{2}(:,1:2);
    bnd_ORTH = EF{tp}{3}(:,1:2);

   
    %%

    [pPD,S_PD] = polyfit(dRES(G),d_PD(G),1); 
    [pOD,S_OD] = polyfit(dRES(G),d_OD(G),1);
    [pORTH,S_ORTH] = polyfit(dRES(G),d_ORTH(G),1);
    
    [qPD,Sq_PD] = polyfit(dslopes(G),d_PD(G),1); 
    [qOD,Sq_OD] = polyfit(dslopes(G),d_OD(G),1);
    [qORTH,Sq_ORTH] = polyfit(dslopes(G),d_ORTH(G),1);

    R2_PD = 1 - sum((polyval(pPD,dRES(G)) - d_PD(G)).^2)./sum((mean(d_PD(G)) - d_PD(G)).^2);
    R2_OD = 1 - sum((polyval(pOD,dRES(G)) - d_OD(G)).^2)./sum((mean(d_OD(G)) - d_OD(G)).^2);
    R2_ORTH = 1 - sum((polyval(pORTH,dRES(G)) - d_ORTH(G)).^2)./sum((mean(d_ORTH(G)) - d_ORTH(G)).^2);

%     fitE_PD = sqrt(diag((S_PD.R)\inv(S_PD.R'))./S_PD.normr.^2./S_PD.df);
%     fitE_ORTH = sqrt(diag((S_ORTH.R)\inv(S_ORTH.R'))./S_ORTH.normr.^2./S_ORTH.df);
%     fitE_OD = sqrt(diag((S_OD.R)\inv(S_OD.R'))./S_OD.normr.^2./S_OD.df);

    end
    
    allcoefs = corrcoef([dRES(G) dslopes(G) d_PD(G) d_ORTH(G) d_OD(G)]);
    
    RESslope_slope = polyfit(dRES(G),dslopes(G),1);
    
    dRESslopes(bsamp,:) = [RESslope_slope(1) pPD(1) pORTH(1) pOD(1)];
    dslopesslopes(bsamp,:) = [RESslope_slope(1) qPD(1) qORTH(1) qOD(1)];

    ccf_dRES(bsamp,:) = allcoefs(1,2:end);
    ccf_dslopes(bsamp,:) = allcoefs(2,3:end);
    
end
%%
[PD1r] = polyfit(ccf_dRES(:,1),ccf_dslopes(:,1),1);
[PD2r] = polyfit(ccf_dRES(:,1),ccf_dRES(:,2),1);
[ORTH1r] = polyfit(ccf_dRES(:,1),ccf_dslopes(:,2),1);
[ORTH2r] = polyfit(ccf_dRES(:,1),ccf_dRES(:,3),1);
[OD1r] = polyfit(ccf_dRES(:,1),ccf_dslopes(:,3),1);
[OD2r] = polyfit(ccf_dRES(:,1),ccf_dRES(:,4),1);
%%
CCFRES = ccf_dRES;
CCFSLP = ccf_dslopes;
xbnd = [-1 0];
ybnd = [-1 1];

figure; hold on; 
subplot(1,3,1); hold on; 
plot(CCFRES(:,1),CCFSLP(:,1),'k.'); plot(CCFRES(:,1),CCFRES(:,2),'r.');
plot(xbnd,polyval(PD1r,xbnd),'k'); plot(xbnd,polyval(PD2r,xbnd),'r');
xlim(xbnd); ylim(ybnd);

subplot(1,3,2); hold on; 
plot(CCFRES(:,1),CCFSLP(:,2),'k.'); plot(CCFRES(:,1),CCFRES(:,3),'r.');
plot(xbnd,polyval(ORTH1r,xbnd),'k'); plot(xbnd,polyval(ORTH2r,xbnd),'r');
xlim(xbnd); ylim(ybnd);

subplot(1,3,3); hold on; 
plot(CCFRES(:,1),CCFSLP(:,3),'k.'); plot(CCFRES(:,1),CCFRES(:,4),'r.');
plot(xbnd,polyval(OD1r,xbnd),'k'); plot(xbnd,polyval(OD2r,xbnd),'r');
xlim(xbnd); ylim(ybnd);

%%
% figure; hold on; 
% subplot(1,3,1); hold on; 
% plot(ccf_dRES(:,1),dRESslopes(:,2),'k.'); plot(ccf_dRES(:,1),dslopesslopes(:,2),'r.');
% % xlim([0 1]); ylim([0 1]);
% 
% subplot(1,3,2); hold on; 
% plot(ccf_dRES(:,1),dRESslopes(:,3),'k.'); plot(ccf_dRES(:,1),dslopesslopes(:,3),'r.');
% % xlim([0 1]); ylim([0 1]);
% 
% subplot(1,3,3); hold on; 
% plot(ccf_dRES(:,1),dRESslopes(:,4),'k.'); plot(ccf_dRES(:,1),dslopesslopes(:,4),'r.');
% % xlim([0 1]); ylim([0 1]);
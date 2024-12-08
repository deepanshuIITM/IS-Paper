function [bPf] = bScarceIS(r, c, nboot,type)
% Generate bootstrap samples for R and C
[~,bootSam_r] = bootstrp(nboot,[],r);
[~,bootSam_c] = bootstrp(nboot,[],c);

bPf = zeros(nboot,1);
for i = 1:nboot
    
    br = r(bootSam_r(:,i));
    bc = c(bootSam_c(:,i));
    
    %%
    % Finding failure Response(R) points
    brcdfC = CDF_TPNT(bc,br);
    brcdfR = CDF_TPNT(br,br);
    
    brG = brcdfC - (1-brcdfR); brI = brG<0;
    
    
    bindexR = find(brI);
    bfailRs = zeros(numel(bindexR),1);
    for bm = 1:numel(bindexR)
        bfailRs(bm,1) = br(bindexR(bm));
    end
    
    % Finding failure Capacity(C) points
    bccdfC = CDF_TPNT(bc,bc);
    bccdfR = CDF_TPNT(br,bc);
    
    bcG = bccdfC - (1-bccdfR); bcI = bcG<0; bcInflex = sum(bcI);
    
    if bcInflex > 0
        bindexC = find(bcI);
        bfailCs = zeros(numel(bindexC),1);
        for bk = 1:numel(bindexC)
            bfailCs(bk,1) = bc(bindexC(bk));
        end
        binta = max(max(bfailRs),max(bfailCs));
        bsafeCs = setdiff(bc,bfailCs);
        bintb = min(bsafeCs);
        binterval = binta:(bintb-binta)/100:bintb;
        %     bmuV = max(max(bfailCs),max(bfailRs));
    else
        binta = max(bfailRs);
        bintb = min(bc);
        binterval = binta:(bintb-binta)/100:bintb;
        %     bmuV = (min(bc)+max(bfailRs))/2;
    end
    binterval = binterval';
    brint = CDF_TPNT(br,binterval);
    bcint = CDF_TPNT(bc,binterval);
    
    bG = bcint-(1-brint); bindX = find(bG<=0);
    bnewFails = zeros(numel(bindX),1);
    for bq = 1:numel(bnewFails)
        bnewFails(bq,1) = binterval(bindX(bq));
    end
    %%
    bmuV = max(bnewFails);
    
    bsdV = (1/2)*(mean(bc)-mean(br));
    
    bv = normrnd(bmuV,bsdV,1e3,1);
    bVpdf = normpdf(bv,bmuV,bsdV);
    %%
    if type == 'int1'
        
        bCcdf = CDF_TPNT(bc,bv);
        bRpdf = akde1d(br,bv);
        bPf(i,:) = mean(bRpdf.*bCcdf./bVpdf);
    else
        bCpdf = akde1d(bc,bv);
        bRcdf = CDF_TPNT(br,bv);
        bPf(i,:) = mean((1-bRcdf).*bCpdf./bVpdf);
    end
end
end
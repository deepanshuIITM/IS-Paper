function [Pf, bootPf, o_count] = scarceIS(r,c,Nis)
%%

nboot = 100;

o_count = 0;

ku = ceil(0.1*Nis);

rcdfC = CDF_TPNT(c,r);
rcdfR = CDF_TPNT(r,r);

rG = rcdfC - (1-rcdfR); rI = rG<0;

% Finding failure R points
indexR = find(rI);
failRs = zeros(numel(indexR),1);
for m = 1:numel(indexR)
    failRs(m,1) = r(indexR(m));
end

ccdfC = CDF_TPNT(c,c);
ccdfR = CDF_TPNT(r,c);

cG = ccdfC - (1-ccdfR); cI = cG<0; cInflex = sum(cI);

% Finding failure C points
if cInflex > 0
    indexC = find(cI);
    failCs = zeros(numel(indexC),1);
    for k = 1:numel(indexC)
        failCs(k,1) = c(indexC(k));
    end
    inta = max(max(failRs),max(failCs));
    safeCs = setdiff(c,failCs);
    intb = min(safeCs);
    interval = inta:(intb-inta)/100:intb;
else
    inta = max(failRs);
    intb = min(c);
    interval = inta:(intb-inta)/100:intb;
end
interval = interval';
rint = CDF_TPNT(r,interval);
cint = CDF_TPNT(c,interval);

G = cint-(1-rint); indX = find(G<=0);
newFails = zeros(numel(indX),1);
for q = 1:numel(newFails)
    newFails(q,1) = interval(indX(q));
end
% MPP of ISD
muV = max(newFails);
% Standard deviation of ISD
sdV = (1/2)*(mean(c)-mean(r));

v = normrnd(muV,sdV,1e3,1);
Vpdf = normpdf(v,muV,sdV);


%% Tail index estimation
%% Upper tail index
xu = sort(r-mean(r),'descend');

Hr = mean(log(xu(1:ku)/xu(ku+1)));

%% Lower tail index
xl = abs(sort(c-mean(c)));

Hc = mean(log(xl(1:ku)/xl(ku+1)));

if Hr>Hc
    o_count = o_count+1;
    Rcdf = CDF_TPNT(r,v);
    Cpdf = akde1d(c,v);
    % Failure probability estimate from original samples
    Pf = mean((1-Rcdf).*Cpdf./Vpdf);
    
    bootstat = bScarceIS(r,c,nboot,'int2');
    % Failure probability estimates from bootstrap samples
    bootPf = bootstat;
else
    
    Ccdf = CDF_TPNT(c,v);
    Rpdf = akde1d(r,v);
    % Failure probability estimate from original samples
    Pf = mean(Rpdf.*Ccdf./Vpdf);
    
    bootstat = bScarceIS(r,c,nboot,'int1');
    % Failure probability estimates from bootstrap samples
    bootPf = bootstat;
end

end


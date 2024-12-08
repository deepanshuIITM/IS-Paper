Nis = 50;
LQ = zeros(3,3);
Qmean = zeros(3,3);
Qsd = zeros(3,3);
%% v critical - 8 (2.94), 16 (3.68), 64 (4.62)
vcr = [4,4.5,5].*1e-3;
% act_beta = [2.98,3.50,3.97];
o_count = zeros(100,3);
for j = 1:3
    Pf = zeros(100,1);
    bootPf = zeros(100,100);
    % Counter set
    fprintf('\nIteration: ')
    for i = 1:100
        %% base random variables
        q = normrnd(20*1e3,2*1e3,Nis,1);
        L = normrnd(6,0.3,Nis,1);
        d = normrnd(25e-2,0.5e-2,Nis,1);
        bf = normrnd(25e-2,0.5e-2,Nis,1);
        tw = normrnd(2e-2,0.2e-2,Nis,1);
        tf = normrnd(2e-2,0.2e-2,Nis,1);
        
        E = normrnd(210e9,10e9,Nis,1);
        
        x = 0.5525*L;
        %% response
        I = (bf.*(d.^3) - (bf-tw).*((d - 2*tf).^3))/12;
        r = (q.*(x.^2).*(4*L.^3 - 8*(L.^2).*x + 5*L.*x.^2 - x.^3))./(120.*L.*I);
        %% capacity
        c = vcr(j)*E;
        
        
        [Pf(i,1), bootPf(i,:), o_count(i,j)] = scarceIS(r,c,Nis);
        % Iteration counter update
        if i>1
            for k=0:log10(i-1)
                fprintf('\b'); % delete previous counter display
            end
        end
        fprintf('%d', i);
    end
    % Percentiles from original samples
    beta = -norminv(Pf);
    Q = quantile(beta, [0.25, 0.5, 0.75]);
%     LQ(j,:) = Q/act_beta(j);
    LQ(j,:) = Q;
    
    % Mean and standard deviation of bootstrap percentiles
    bbeta = -norminv(bootPf);
    Qboot = quantile(bbeta, [0.25, 0.5, 0.75], 2);
%     LQboot = Qboot/act_beta(j);
    LQboot = Qboot;
    
    Qmean(j,:) = mean(LQboot,1);
    Qsd(j,:) = std(LQboot,0,1);
end
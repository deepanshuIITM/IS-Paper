% N = 1e7;
% u1 = normrnd(0,1,N,1);
% u2 = normrnd(0,1,N,1);
% 
% G = -2.62 + (0.15*u1.^2) + u2;
% 
% % pf = 0.9916
% pf = mean(G<=0);
% 
% norminv(pf)
Nis = 200;
LQ = zeros(3,1);
Qmean = zeros(3,1);
Qsd = zeros(3,1);

o_count = zeros(100,1);
for j =1:1
    Pf = zeros(100,1);
    bootPf = zeros(100,100);
    % Counter set
    fprintf('\nIteration: ')
    for i = 1:100
        %% response
        r =  0.15*normrnd(0,1,Nis,1).^2;
        %% capacity
        c = 2.62 - normrnd(0,1,Nis,1);
        
        
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
    LQ = Q;
    
    % Mean and standard deviation of bootstrap percentiles
%     bbeta = -norminv(bootPf);
    Qboot = quantile(bootPf, [0.25, 0.5, 0.75], 2);
%     LQboot = Qboot/act_beta(j);
    LQboot = -norminv(Qboot);
    
    Qmean = mean(LQboot,1);
    Qsd = std(LQboot,0,1);
end
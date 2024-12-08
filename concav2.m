% N = 1e7;
% x1 = normrnd(0,1,N,1);
% x2 = normrnd(0,1,N,1);
%
% G = -5*x1.^2 + x2.^2 +45;
%
% % pf = 0.0024
% pf = mean(G<=0);
%
% norminv(pf)


Nis = 50;

o_count = zeros(100,1);
t_beta = norminv(1-0.0024);

Pf = zeros(100,1);
bootPf = zeros(100,100);
% Counter set
fprintf('\nIteration: ')
for i = 1:100
    %% response
    r =  5*normrnd(0,1,Nis,1).^2;
    %% capacity
    c = normrnd(0,1,Nis,1).^2 + 45;
    
    
    [Pf(i,1), bootPf(i,:), o_count(i,1)] = scarceIS(r,c,Nis);
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
%   LQ(j,:) = Q/act_beta(j);
LQ = Q/t_beta;

% Mean and standard deviation of bootstrap percentiles
bbeta = -norminv(bootPf);
Qboot = quantile(bbeta, [0.25, 0.5, 0.75], 2);
%     LQboot = Qboot/act_beta(j);
%     LQboot = Qboot;

Qmean = mean(Qboot,1);
Qs = std(Qboot,0,1);
LQmean = Qmean./t_beta;
LQsd = Qs./t_beta;

% N = 1e9;
% q = normrnd(20000,1600,N,1);
% l = normrnd(12,0.24,N,1);
% As = normrnd(9.82e-4,5.89e-5,N,1);
% Ac = normrnd(0.04,0.008,N,1);
% Es = normrnd(1.2e11,8.4e9,N,1);
% Ec = normrnd(3e10,2.4e9,N,1);
%
% Z = 0.03;
%
% G = (Z./q)-((0.5*l.^2).*((3.81./(Ac.*Ec)) + (1.13./(As.*Es))));
% G = (Z)-((0.5*q.*l.^2).*((3.81./(Ac.*Ec)) + (1.13./(As.*Es))));
%
% % pf =5.0215e-4 beta=3.2893 with MCS sample size = 1e9
% pf = mean(G<0);
%
% norminv(1-pf)


Nis = 50;

o_count = zeros(100,1);
t_beta = norminv(1-5.0215e-4);

Pf = zeros(100,1);
bootPf = zeros(100,100);
% Counter set
fprintf('\nIteration: ')
for i = 1:100
    q = normrnd(20000,1600,Nis,1);
    l = normrnd(12,0.24,Nis,1);
    As = normrnd(9.82e-4,5.89e-5,Nis,1);
    Ac = normrnd(0.04,0.008,Nis,1);
    Es = normrnd(1.2e11,8.4e9,Nis,1);
    Ec = normrnd(3e10,2.4e9,Nis,1);
    
    Z = 0.03;
    %% response
    r = ((0.5*l.^2).*((3.81./(Ac.*Ec)) + (1.13./(As.*Es))));
    %% capacity
    c = (Z./q);
    
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
%     LQ(j,:) = Q/act_beta(j);
LQ = Q/t_beta;

% Mean and standard deviation of bootstrap percentiles
bbeta = -norminv(bootPf);
Qboot = quantile(bbeta, [0.25, 0.5, 0.75], 2);
%     LQboot = Qboot/act_beta(j);
LQboot = Qboot;

Qmean = mean(LQboot,1);
Qsd = std(LQboot,0,1);

LQmean = Qmean/t_beta;
LQsd = Qsd/t_beta;

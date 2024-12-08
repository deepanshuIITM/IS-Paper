Nis = 1e2;
%
% q = normrnd(20000,1600,Nis,1);
% l = normrnd(12,0.24,Nis,1);
% As = normrnd(9.82e-4,5.89e-5,Nis,1);
% Ac = normrnd(0.04,0.008,Nis,1);
% Es = normrnd(1.2e11,8.4e9,Nis,1);
% Ec = normrnd(3e10,2.4e9,Nis,1);
%
% Z = 0.03;
% %% response
% r = ((0.5*l.^2).*((3.81./(Ac.*Ec)) + (1.13./(As.*Es))));
% %% capacity
% c = (Z./q);

%% response
r =  5*normrnd(0,1,Nis,1).^2;
%% capacity
c = normrnd(0,1,Nis,1).^2 + 45;

figure()
plot(r,ksdensity(r,r),'r.','MarkerSize',7)
hold on
plot(c,ksdensity(c,c),'b.','MarkerSize',7)
hold off
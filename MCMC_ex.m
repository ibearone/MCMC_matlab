clc;clear;close all;

%%% Define Parameter
L=2; J=1; MC_circle=1000;h=1;
Tmax=5;
T_set=6
T=(1:1:T_set)/T_set*5;
for k=1:T_set
Gamma{1}{k}=randi([0 1], L)*2-1;
end
for n=2:MC_circle
[Gamma{n}] = MCMC_metropolis_single(Gamma{n-1},L,T,J,h,1,mod(n,2))
end
%%%reshape
for k=1:T_set
    for n=1:MC_circle
        Gamma_rs{n,k}=Gamma{n}{k}
        [E_tot(n,k),Ms_tot(n,k)] = E_Ms_int_single(Gamma_rs{n,k},L,T(k),J,h)
            E_tot_2(n,k)=1/n*sum(E_tot(1:n,k))/L^2
            Ms_tot_2(n,k)=1/n*sum(Ms_tot(1:n,k))/L^2
    end
end

for k=1:T_set
plot(Ms_tot_2(:,k))
hold on;
end
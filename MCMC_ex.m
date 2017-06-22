%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example simulation with exchange replica
% 2017/6/23
% Guanxiong Qu
% quguanxiong@gmail.com
%
% Parameters:
% Gamma_old: old configuration
% L:         Dimension
% T:         Temperature
% J:         Exchange
% h:         bias field
% ex:        exchange replica (0=fasle/ 1=ture)
% eo:        exchange sequence (0=even/ 1=odd)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

%%% Define Parameter
L=2; J=1; MC_circle=1000;h=1;ex=1;Tmax=5;T_set=6
%%% Initializing T set
T=(1:1:T_set)/T_set*Tmax;
%%% Generating corresponding configuration set
for k=1:T_set
Gamma{1}{k}=randi([0 1], L)*2-1;
end
%%% MCMC with exchange replica
for n=2:MC_circle
[Gamma{n}] = MCMC_metropolis_single(Gamma{n-1},L,T,J,h,ex,mod(n,2))
end
%%% reshape&properties calcuations
for k=1:T_set
    for n=1:MC_circle
        Gamma_rs{n,k}=Gamma{n}{k}
        [E(n,k),Ms(n,k)] = E_Ms_int_single(Gamma_rs{n,k},L,T(k),J,h)
            E_tot(n,k)=1/n*sum(E(1:n,k))/L^2
            Ms_tot(n,k)=1/n*sum(Ms(1:n,k))/L^2
    end
end
%%%plot
subplot(1,2,1)
for k=1:T_set
plot(Ms_tot(:,k))
hold on;
end
%%%plot
subplot(1,2,2)
for k=1:T_set
plot(E_tot(:,k))
hold on;
end

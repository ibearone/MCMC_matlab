%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single run MCMC_metropolis(exchange replica)
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
function[Gamma_new,r,P_ratio] = MCMC_metropolis_single(Gamma_old,L,T,J,h,ex,eo)

%%% exchange replica
if ex==0 %false
%%%initialize
Index_ij=randi([1 L], 1,2);
r=rand(1,1);
%%% dummy upate
Gamma_new=Gamma_old
    i=Index_ij(1);j=Index_ij(2)
    Gamma_new(i,j)=-1*Gamma_new(i,j)
%%% Energy Difference
    [E_new] = E_Ms_int_single(Gamma_new,L,T,J,h)
    [E_old] = E_Ms_int_single(Gamma_old,L,T,J,h)
    Delta_E=E_new-E_old
%%% Probability
    P_ratio=exp(-Delta_E/T)
%%% Metropolis criterion 
    if r<=P_ratio
    Gamma_new=Gamma_new
    else
    Gamma_new=Gamma_old
    end
    
    %%%exchange replica
else  %%%ture
%%%initialize
    rep_num=length(T)
    Index_ij=randi([1 L],rep_num,2);
    r=rand(rep_num,1);
    Gamma_new=Gamma_old
%%% dummy update
    for k=1:rep_num
    Gamma_new{k}(Index_ij(k,1),Index_ij(k,2))=-1*Gamma_new{k}(Index_ij(k,1),Index_ij(k,2))
%%% Energy Difference
    [E_old(k),Ms] = E_Ms_int_single(Gamma_old{k},L,T(k),J,h)
    [E_new(k),Ms] = E_Ms_int_single(Gamma_new{k},L,T(k),J,h)
    Delta_E(k)=E_new(k)-E_old(k)
%%% Probability
    P_ratio(k)=exp(-Delta_E(k)/T(k))
%%% Metropolis criterion 
    if r(k)<=P_ratio(k)
    Gamma_new{k}=Gamma_new{k}
    else
    Gamma_new{k}=Gamma_old{k}
    end
    end
%%% exchnage among T
    r_ex=rand(rep_num/2,1);
%%% even exchange 
    if eo==0
    for n=1:rep_num/2
%%% Metropolis criterion 
    P_ratio_ex(n)=exp(-(E_new(2*n)-E_new(2*n-1))./(T(2*n)-T(2*n-1)))
    if r_ex(n)<=P_ratio_ex(n)
    A=Gamma_new{2*n}
    B=Gamma_new{2*n-1}
    Gamma_new{2*n}=B
    Gamma_new{2*n-1}=A
    else
    end 
    end
    else
%%% odd exchange
    for n=1:rep_num/2-1
    P_ratio_ex(n)=exp(-(E_new(2*n+1)-E_new(2*n))/(T(2*n+1)-T(2*n)))
    if r_ex(n)<=P_ratio_ex(n)
    A=Gamma_new{2*n+1}
    B=Gamma_new{2*n}
    Gamma_new{2*n+1}=B
    Gamma_new{2*n}=A
    else
    end 
    end
    end
    end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single run Energy and Magnetization
% 2017/6/23
% Guanxiong Qu
% quguanxiong@gmail.opcm
%
% Parameters:
% Gamma: configuration
% L:         Dimension
% T:         Temperature
% J:         Exchange
% h:         bias field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E,Ms] = E_Ms_int_single(Gamma,L,T,J,h)
%%%Energy 
    ij=0
    for i=1:L
        for j=1:L-1
        ij=ij+Gamma(i,j)*Gamma(i,j+1)
        end
    end
        for j=1:L
        for i=1:L-1
        ij=ij+Gamma(i,j)*Gamma(i+1,j)
        end
    end
    E=-h*sum(sum(Gamma))-J*ij
%%% Magnetization
    Ms=sum(sum(Gamma))
end

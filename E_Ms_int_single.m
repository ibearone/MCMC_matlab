function [E_2,Ms] = E_Ms_int_single(Gamma,L,T,J,h)
    %%%E
    E=0
    for i=1:L
        for j=1:L-1
        E=E+Gamma(i,j)*Gamma(i,j+1)
        end
    end
        for j=1:L
        for i=1:L-1
        E=E+Gamma(i,j)*Gamma(i+1,j)
        end
    end
    E_2=-h*sum(sum(Gamma))-J*E

    Ms=sum(sum(Gamma))
end
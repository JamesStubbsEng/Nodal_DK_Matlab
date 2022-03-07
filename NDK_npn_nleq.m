function [in, jacobian] = bjt_nleq(v)

%% BJT Parameters:
Is = 1.16e-14;
Bf = 200;
Br = 3;
Vt = 0.02528;

%BJT currents
in = Is*[   
        1/Bf * (exp((v(2)-v(1))/Vt) - 1) + 1/Br * (exp(-v(1)/Vt) - 1);
        -(-(exp(-v(1)/Vt)-1) + (1+Bf)/Bf * (exp((v(2)-v(1))/Vt) - 1));];

jacobian = Is/Vt*[ 
        -1/Bf * exp((v(2)-v(1))/Vt) - 1/Br * exp(-v(1)/Vt),     1/Bf * exp((v(2)-v(1))/Vt);
        -exp(-v(1)/Vt) + (1 + Bf)/Bf * exp((v(2)-v(1))/Vt),     -(1 + Bf)/Bf * exp((v(2)-v(1))/Vt)];
end


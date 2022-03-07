function [in, jacobian] = diode_nleq(vn)
%Non-linear elements parameters

% Be careful to make the NL currents in the same direction as DK
% definitions!!!! In this case, the diode current is made 
% negative.

% diode parameters
Is = 1e-15;
Vt = 25.85e-3;
eta = 1;
%diode current
in = -2*Is * sinh(vn/(eta*Vt));
jacobian = -2*Is/(eta*Vt) * cosh(vn/(eta*Vt));

end


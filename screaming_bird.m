%simulation parameters
Fs = 44100.0;
T = 1/Fs;

% needs manual input
% -------------------------------------------
Vcc = 9;

%Resistors values
R1 = 1e6;
R2 = 43e3;
R3 = 430e3;
R4 = 390;
R5 = 10e3;
P1 = 100e3;
Rl = 200e3;
%Capacitors and Inductors values
C1 = 2.2e-9;
C3 = 2.2e-9;
%Variable Resistors values
NR = sparse([...
% 1 2 3 4 5 6
1, 0, 0, 0, 0, 0;...% R1
0, 1, 0, 0, 0, 0;...% R2
0, -1, 1, 0, 0, 0;...% R3
0, 0, 0, 0, 1, 0;...% R4
0, 0, 1, -1, 0, 0;...% R5
0, 0, 0, 0, 0, 1;...% P1
0, 0, 0, 0, 0, 1 ...% Rl
]);

Nx = sparse([...
% 1 2 3 4 5 6
1, -1, 0, 0, 0, 0;...% C1
0, 0, 0, 1, 0, -1 ...% C3
]);
Nv = sparse([...
]);
Nn = sparse([...
% 1 2 3 4 5 6
0, 1, 0, 0, -1, 0;...% T1
0, 0, 0, 1, -1, 0 ...% T1
]);
Nu = sparse([...
% 1 2 3 4 5 6
1, 0, 0, 0, 0, 0;...% Vin
0, 0, 1, 0, 0, 0 ...% Vcc
]);
NuT = sparse([...
% 1 2
1, 0;...
0, 0;...
0, 1;...
0, 0;...
0, 0;...
0, 0 ...
]);
No = sparse([...
% 1 2 3 4 5 6
0, 0, 0, 0, 0, 1 ...% Vout
]);
GR = sparse(diag([1/R1,1/R2,1/R3,1/R4,1/R5,1/P1,1/Rl]));
Gx = sparse(diag([2*C1/T,2*C3/T]));
Z = sparse(diag([ 1, 1]));
Vr = [];
vnin = [];
[rowsNn colsNn] = size(Nn);
Vn = zeros(rowsNn,1);

[rowsNx colsNx] = size(Nx);
x = zeros(rowsNx,1);

% -------------------------------------------

% compute system matrix S

S11 = NR'*GR*NR + Nx'*Gx*Nx
S12 = Nu'
S21 = Nu
S22 = zeros(size(Nu,1))

S = [S11 S12; S21 S22]

Si = inv(S);

% padded matrices
Nrp = [NR zeros(size(NR,1),size(Nu,1))]
Nxp = [Nx zeros(size(Nx,1),size(Nu,1))]
Nnp = [Nn zeros(size(Nn,1),size(Nu,1))]
Nop = [No zeros(size(No,1),size(Nu,1))]
Nup = [Nu zeros(size(Nu,1),size(Nu,1))]  
% resolve voltage sources matrix Nup to proper format under new variable Nup2  
Nup2 = [zeros(size(NR,2),size(Nu,1)); eye(size(Nu,1)); zeros(0,size(Nu,1))]

% compute state-space system matrices

A = 2*Z*Gx*Nxp*Si*Nxp' - Z
B = 2*Z*Gx*Nxp*Si*Nup2
C = 2*Z*Gx*Nxp*Si*Nnp'

D = Nop*Si*Nxp'
E = Nop*Si*Nup2
F = Nop*Si*Nnp'

G = Nnp*Si*Nxp'
H = Nnp*Si*Nup2
K = Nnp*Si*Nnp'

% Setup input signal
t = [0:T:0.1].';
Vi = 1*sin(2*pi*400*t);

% Initialize Output Signal
N = length(Vi);
Vo = zeros(N,1);


for n = 1:N
    u = [Vi(n,1) Vcc]';

    p = G*x + H*u;

    [Vn, In] = nonlinearSolver(p,K,Vn, @sb_nleq);

    Vo(n,1) = D*x + E*u + F*In;
    x = A*x + B*u + C*In;

end

plot(t,Vi,t,Vo);
axis([0 0.1 -9 9]);


%simulation parameters
Fs = 44100.0;
T = 1/Fs;

% needs manual input
% -------------------------------------------
% passive components
R1 = 240e3;
C1 = 80e-12;

% diode parameters
Is = 1e-15;
Vt = 25.85e-3;
eta = 1;


% connection matrices
NR = [1 -1];
Nv = [0];
Nx = [0 1];
Nu = [1 0];
Nn = [0 1];

GR = diag([1/R1]);
Rv = [0];
Gx = diag([2*C1/T]);

No = [0 1];

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
% needs manual input
% -------------------------------------------
Z = diag([1]);
% -------------------------------------------

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

% Initialize State and intermediate vectors

[rowsNn colsNn] = size(Nn);
Vn = zeros(rowsNn,1);

[rowsNx colsNx] = size(Nx);
x = zeros(rowsNx,1);

for n = 1:N
    u = Vi(n,1);

    p = G*x + H*u;

    [Vn, In] = nonlinearSolver(p,K,Vn, @diode_nleq);

    Vo(n,1) = D*x + E*u + F*In;
    x = A*x + B*u + C*In;

end

plot(t,Vi,t,Vo);
axis([0 0.1 -9 9]);


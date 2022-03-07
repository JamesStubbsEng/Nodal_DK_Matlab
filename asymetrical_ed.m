%simulation parameters
Fs = 100000.0;
T = 1/Fs;

% needs manual input
% -------------------------------------------
Vcc = 9;

% passive components
R1 = 100e3;
R2 = 22;
R3 = 470e3;
R4 = 10e3;
R5 = 700;
C1 = 0.047e-6;
C2 = 0.25e-9;
C3 = 0.47e-6;


% connection matrices
NR = [0 0 0 0 1 0; 0 0 1 0 0 0; 0 0 0 0 -1 1; 0 1 0 0 0 -1; 0 0 0 1 0 0];
Nv = [0];
Nx = [-1 0 0 0 1 0; 0 0 0 0 1 -1; 0 0 0 -1 0 1];
Nu = [1 0 0 0 0 0; 0 1 0 0 0 0];
Nn = [0 0 0 0 -1 1; 0 0 -1 0 0 1];

GR = diag([1/R1 1/R2 1/R3 1/R4 1/R5]);
Rv = [0];
Gx = diag([2*C1/T 2*C2/T 2*C3/T]);

No = [0 0 0 1 0 0];

% ------------------------------------------

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

% Needs manual input

% ----------------------------------------
Z = diag([1 1 1]);
% ----------------------------------------

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
Vi = 0.95*sin(2*pi*1000*t);

% Initialize Output Signal
N = length(Vi);
Vo = zeros(N,1);

% Initialize State and intermediate vectors
[rowsNn colsNn] = size(Nn);
Vn = zeros(rowsNn,1);

[rowsNx colsNx] = size(Nx);
x = zeros(rowsNx,1);

for n = 1:N
    u = [Vi(n,1) Vcc]';

    p = G*x + H*u;

    [Vn, In] = nonlinearSolver(p,K,Vn, @NDK_npn_nleq);

    Vo(n,1) = D*x + E*u + F*In;
    x = A*x + B*u + C*In;

end

plot(t,Vi,t,Vo);
axis([0 0.1 -9 9]);



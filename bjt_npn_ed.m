%simulation parameters
Fs = 44100.0;
T = 1/Fs;

% needs manual input
% -------------------------------------------

%rail voltage
Vcc = 9;

% passive components
R1 = 470e3;
R2 = 2.2e3;
R3 = 220;
R4 = 50e3;
C1 = 2.2e-9;
C2 = 2.2e-9;




% connection matrices
NR = [0 -1 0 0 0 1; 0 0 -1 0 0 1; 0 0 0 1 0 0; 0 0 0 0 1 0];
Nv = [0];
Nx = [1 -1 0 0 0 0; 0 0 1 0 -1 0];
Nu = [1 0 0 0 0 0; 0 0 0 0 0 1];
Nn = [0 1 0 -1 0 0; 0 0 1 -1 0 0];

GR = diag([1/R1 1/R2 1/R3 1/R4]);
Rv = [0];
Gx = diag([2*C1/T 2*C2/T]);

No = [0 0 0 0 1 0];

% -------------------------------------------

% compute system matrix S

S11 = NR'*GR*NR + Nx'*Gx*Nx;
S12 = Nu';
S21 = Nu;
S22 = zeros(2);

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
Z = diag([1 1]);
% -------------------------------------------
A =  2*Z*Gx*Nxp*Si*Nxp' - Z
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
Vi = 0.9*sin(2*pi*5000*t);

% Initialize Output Signal
N = length(Vi);
Vo = zeros(N,1);

% Initialize State and intermediate vectors
x = [0 0]';
In = [0 0]';
Vn = [0 0]';


for n = 1:N
    u = [Vi(n,1) Vcc]';
    
    p = G*x + H*u;

    [Vn, In] = nonlinearSolver(p,K,Vn,@bjt_nleq);

    Vo(n,1) = D*x + E*u + F*In;
    x = A*x + B*u + C*In;

end

plot(t,Vi,t,Vo);
axis([0 0.001 -1 1]);


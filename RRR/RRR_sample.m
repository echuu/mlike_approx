
%% Reduced rank regression 
%% Y= XAB' + eps; 
%% Y is n*q, X is n*p, A is p*r, B is q*r, eps ~ N(0,sig2*I_q \kron I_n)
%% Prior for B': N(0,sig2*I_q \kron del^{-2}I_r)
%% Prior for A : N(0,sig2*I_r \kron del^{-2}I_p)

% cd("C:\\Users\\ericc\\mlike_approx\\RRR")

pkg load statistics

% Likelihood parameters

p = 2;           % number of columns in X
q = 3;           % number of columns in Y
r = 2;           % number of columns in B and A
% n = 100;         % number of rows in X and Y
% n = p;
sig2 = 10^(-2);  % fixed for now.

% Generate covariates
sigX = eye(p);

X = mvnrnd(zeros(1, p), sigX, n);
X = eye(p);

% csvwrite('X.csv', X)

% Generate and fix parameters A & B
A = normrnd(0, 1, p, r);
B = normrnd(0, 1, q, r); 

% csvwrite('A.csv', A);
% csvwrite('B.csv', B);

% Generate data
eps = normrnd(0, sqrt(sig2), n, q);
Y = X * A * B' + eps;
C = A * B';

% csvwrite('eps.csv', eps)

% Prior parameters
del=10^(-2);

%% Gibbs sampling

nMCMC = 1500;
B = normrnd(0, 1, q, r);

% csvwrite('B_init.csv', B);

D = r * p + q * r;      % dimension of each MC sample, u
u_df = zeros(nMCMC, D); % store each MCMC sample, u,  row-wise in u_df

for g = 1:nMCMC
    
    % sample from Bt | -
    Btvarpart = sig2 * inv(A' * (X' * X) * A + del^2);
    Btvar = kron(eye(q), Btvarpart);
    Btmu = Btvarpart * A' * X' * Y ./ sig2; 
    
    
    Bt_row = mvnrnd(reshape(Btmu, 1, []), Btvar, 1);  % (1 x rq) row vector
    Bt = reshape(Bt_row, r, q);                       % (r x q)  matrix
    B = Bt';
    
    Mt = inv(B' * B) * (B' * Y' * X) * inv(X' * X);
    M = Mt';
    
    BtB = B' * B; 
    XtX = X' * X;
    
    % sample from A | -
    Avarpart = kron(BtB, XtX / sig2);
    Avar = inv(del^2 / sig2 + Avarpart);
    Amu = Avar * Avarpart * reshape(M, 1, [])';
    
    A_row = mvnrnd(Amu, Avar, 1); % (1 x pr) row vector
    A = reshape(A_row, p, r);     % (p x r)  matrix

    % store the row vectors A (1 x pr), B (1 x rq) together to form the sample
    % u (1 x  (pr + rq)) row vector --> this will follow the data structure of
    % u_df as is the resulting data frame in previous sampling schemes

    u_df(g,:) = [A_row, Bt_row];

end

% csvwrite('u_df_rrr.csv', u_df)

A_g  = reshape(u_df(g,1:p*r), p, r);   % recover the matrix A
Bt_g = reshape(u_df(g,p*r+1:D), r, q); % recover the matrix B^T =: B


C

A * B'

% u_df((g-6):g,:)
% u_df(1:6,:)

% writematrix(u_df, 'u_df_rrr.csv')


% A * B'







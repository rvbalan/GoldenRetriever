function [myerr,x,A,y,z,Relerrs,Objective] = WF_fun(x,A,y)
%x is the true signal, n x 1 vector
%A is the frame vectors in rows (m x n matrix), so F' in the GR code
%y (m x 1 observation vector) is optional, if no y is provided then it is computed without noise
%% Make signal and data 

%x = randn(n,1) + 1i*randn(n,1);
if ~exist('x','var')
  n = 5;
  x = (1/sqrt(n))*randn(n,1);
end
n=length(x);
                    
%A = 1/sqrt(2)*randn(m,n) + 1i/sqrt(2)*randn(m,n);
if ~exist('A','var')
  m = round(3*n); 
  A = randn(n,m);
end
A=A';
m=size(A,1);

if ~exist('y','var')
    y = abs(A*x).^2 ;
end

%% Initialization

npower_iter = 50;                           % Number of power iterations 
z0 = randn(n,1); z0 = z0/norm(z0,'fro');    % Initial guess 
for tt = 1:npower_iter,                     % Power iterations
    z0 = A'*(y.* (A*z0)); z0 = z0/norm(z0,'fro');
end

normest = sqrt(sum(y)/numel(y));    % Estimate norm to scale eigenvector  
z = normest * z0;                   % Apply scaling 
Relerrs = norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro'); % Initial rel. error
Objective=norm(y-abs(A*z).^2)^2/m;
%% Loop

T = 25000;                           % Max number of iterations
tau0 = 330;                         % Time constant for step size
mu = @(t) min(1-exp(-t/tau0), 0.002); % Schedule for step size

for t = 1:T,
    yz = A*z;
    grad  = 1/m* A'*( ( abs(yz).^2-y ) .* yz ); % Wirtinger gradient
    z = z - mu(t)/normest^2 * grad;             % Gradient update 

    %Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro')];
    %Objective = [Objective, norm(y-abs(A*z).^2)^2/m];
end
 
%% Check results

%  fprintf('Relative error after initialization: %f\n', Relerrs(1))
%  fprintf('Relative error after %d iterations: %f\n', T, Relerrs(T+1))
 %myerror=Objective(T+1)
 myerr=min(norm(z-x),norm(z+x));
 
%  figure(1), semilogy(0:T,Relerrs) 
%  xlabel('Iteration'), ylabel('Relative error (log10)'), ...
%      title('Relative error vs. iteration count')
%  
%  figure(2), semilogy(0:T,Objective) 
%  xlabel('Iteration'), ylabel('Least Square Error (log10)'), ...
%      title('Least Square Error vs. iteration count')

    
  

end

%Code from Mahdi Soltanolkotabi's website 1D Gaussian


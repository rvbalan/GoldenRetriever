function [x,lambda,counter] = fixedPointAny(xest,lambdaest,param,F,Rzero,c)
%This is a fixed point correction function for marching in the any
%direction.
%xiest is the algorithms estimate to xi
%lambdaest is the algorithms estimate to lambda
%param is a structure which contains, among other things, the sizes m and
%n, as well as tolerance (which says how close successive iterations should
%be) and counter (max number of iterations)
%param also contains the semidefinite matrix Q.
%Czero is a matrix (2nx2n) which is defined in runFunction()
%c is an index 1<=c<=2n+1, which column was switched when this was called,
%defined in stepCFun()

n=param.n;

% if ~exist('c','var')
%     c=n+1;
% else


Q=param.Q;

%Given x,lambda,c, compute x'' (xdp),x' (xp),xc
x=xest;
lambda=lambdaest;
cvect=zeros(1,n+1);
cvect(c)=1;

zeta=[x;lambda];
err=inf;
counter=0;

%Start the loop here
while err>param.tolerance && counter<param.counter
    counter=counter+1;
    %xdp=[zeta(1:c); 0 ; zeta(c+1:end-1)];
    %This creates the truncated normal hessian with appended last column of
    %Qx
    x=zeta(1:n);
    lambda=zeta(n+1);
    myR=Rfun(F,x,param);
    myHess=3*myR + lambda*Q - Rzero;
    A=[myHess, Q*x; cvect]; %Extended Hessian, with last vector as all zeros except for 1
    gradJ=(myR + lambda*Q - Rzero)*x; 
    b=[-gradJ; 0];
    %delta=(A'*A)\(A'*b);
    v=A'*b;
    w=v;
    w(c)=0;
    delta=w/norm(w);
    alpha=Rfun(F,delta(1:n),param)*x;
    beta=myHess*delta(1:n);
    %norm0=norm(gradJ);
    gradJn2=(Rfun(F,x-2*delta(1:n),param) + (lambda-2*delta(n+1))*Q - Rzero)*(x-2*delta(1:n));
    gradJ2=(Rfun(F,x+2*delta(1:n),param) + (lambda+2*delta(n+1))*Q - Rzero)*(x+2*delta(1:n));
    normn2=norm(gradJn2);
    norm2=norm(gradJ2);
    Rdd=Rfun(F,delta(1:n),param)*delta(1:n);
    [step,normgradStep]=PolyMin(alpha,beta,gradJ,Rdd,normn2,norm2,x,lambda,param,Rzero,delta,F);
    %step=directMin(alpha,beta,F,norm0,norm2,x,lambda,param,Rzero,delta)
    zeta2 = zeta+step*delta; % step from line search (of 1/2 * norm(gradJ)^2) with contraint 0<step<2
    %err=norm(zeta2-zeta);
    err=normgradStep;
    zeta=zeta2;
end
%counter
x=zeta(1:n);
lambda=zeta(n+1);
% if counter<10000
%     disp('test here');
% end
%errorX=norm(x-xiest);
%errorLambda=norm(lambda-lambdaest);
end
function [xest,lambdaest,myCount] = fixedPointNewton(xest,lambdaest,param,F,Rzero)
%Does a Newton Corrector as a fixed point iteration

err=inf;
counter=0;
Q=param.Q;
x=xest;
lambda = lambdaest;

while err>param.tolerance && counter<param.counter
    counter=counter+1;
    Rx = Rfun(F,x,param);
    temp = Rx+lambda*Q - Rzero;
    myGradient = temp*x;
    %myFullGrad = [myGradient; (1/2)*x'*Q*x];
    myHessian = 2*Rx + temp;
    fullHess = [myHessian, Q*x];
    errvec=((fullHess)\myGradient);
    x = x - errvec(1:end-1); %Update the x terms
    lambda = lambda - errvec(end); %Update the lambda terms
    err=norm(errvec); %Compute the error of the system
end
xest = x;
lambdaest = lambda; 
myCount= counter;

end


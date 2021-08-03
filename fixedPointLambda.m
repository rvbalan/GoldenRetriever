function [myEstimate,myCount] = fixedPointLambda(xest,lambdaest,param,F,Rzero)
%This is a fixed point correction function for marching in the lambda
%direction.
%xest is the algorithms estimate to x
%lambdaest is the algorithms estimate to lambda
%param is a structure which contains, among other things, the sizes m and
%n, as well as tolerance (which says how close successive iterations should
%be) and counter (max number of iterations)
err=inf;
    counter=0;
    zeta=xest;
    while err>param.tolerance && counter<param.counter
        counter=counter+1;
        myR=Rfun(F,zeta,param); %R(zeta)
        zeta2 = (3*myR+lambdaest*param.Q-Rzero)\((2*myR)*zeta);
        err=norm(zeta2-zeta);
        zeta=zeta2;
    end
     myEstimate=zeta;
     myCount=counter;
end



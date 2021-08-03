function [myNorm,myIX,mygradJ,myJX]=calcDebuger1(xest,y,F,param,Rzero,lambdaest)

mygradJ=norm((Rfun(F,xest,param) + lambdaest*param.Q -Rzero)*xest); %Norm of the gradient
%mygradI=norm((Rfun(F,xest,param) -Rzero)*xest); %Norm of the gradient
myNorm=norm(xest)^2; %Norm of x^2
mySum=0;
for k=1:param.m
    mySum=mySum+(y(k)-abs(dot(xest,F(:,k)))^2)^2;
end
myIX=(1/(4*param.m))*mySum; %Ix Criterion

myJX=myIX+lambdaest/2 * dot(param.Q*xest,xest); %Jx Criterion
end


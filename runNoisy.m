function [myEstimate] = runNoisy(y,F,z)
%RUNNOISY This function attempts to find the global minimum to the noisy quadratic
%loss, it does so using a homotopy algorithm from the global minimum to the
%non-noisy quadratic loss
%z is the clean signal, nx1
%F is the frame vectors, arranged in an nxm matrix
%y is the noisy measumerments, mx1


%Initialize Parameters
n=length(z);
m=length(y);
h=1e-3; %step size
totsteps=10^5; %Total number of steps in homotopy path
mythreshold=1e-9; %Threshold for fixed point correction
switchDirectionThreshold=100; %threshold for condition number to switch to another direction

%Create R(z)
Rz=Rfun(F,z);
%Create Rzero
Z = F.*(ones(n,1)*(y'));
Rzero = (Z*F')/m;


%Initialize the lambda and x parameters
lambda=1; %lambda goes from lambda=1 to (hopefully) lambda=0
xest=z; %initialized at x=z, the clean signal

%The extended estimate vector [x,lambda]
myEst=[xest;lambda];

%Now we Initialize the directions, marching down along lambda
dx=2*Rz\((Rz-Rzero)*z);
dlambda=-1;
dEst=[dx;dlambda];

%Until lambda=0 or total steps is exceeded, we march along the path
numsteps=0;
while lambda>0 && numsteps<totsteps
    %Increment steps
    numsteps=numsteps+1;
    
    %Test if you are stepping past lambda=0
    if myEst(n+1)+h*dEst(n+1)<0
        %If past lambda=0, step size adjusted to get to exactly lambda=0
        hp=-myEst(n+1)/dEst(n+1);
        myEst=myEst+hp*dEst;
        break
    end
    
    tempEst=myEst+h*dEst; %We do a step along the tangent direction
    mysign=sign(tempEst-myEst);%See how the signs of the components changed
    myEst=tempEst;
    
    %Now we do a fixed point correction
    error=inf;
    numfixedsteps=0;
    while error>mythreshold && numfixedsteps<totsteps/10^3
        numfixedsteps=numfixedsteps+1;
        x=myEst(1:n);
        lambda=myEst(n+1);
        Rx=Rfun(F,x);
        %Hessext=[3*Rx-lambda*Rz-(1-lambda)*Rzero,-(Rz-Rzero)*x;[zeros(1,n),1]];
        Hessext=[3*Rx-lambda*Rz-(1-lambda)*Rzero,-(Rz-Rzero)*x];
        gradF=Rx*x-((lambda*Rz)+(1-lambda)*Rzero)*x;
        error=norm(gradF);
        %myEst=myEst-Hessext\[gradF;0];
        myEst=myEst-Hessext\gradF;
    end
    x=myEst(1:n);
    lambda=myEst(n+1);
    Rx=Rfun(F,x);
    Hessext=[3*Rx-lambda*Rz-(1-lambda)*Rzero,-(Rz-Rzero)*x];
    
    %Now you compute the direction of the tangent from the extended hessian
    mydir=null(Hessext);
    %Normalize it to be +/- 1 along the c component 
    if cond(3*Rx-lambda*Rz-(1-lambda)*Rzero)>switchDirectionThreshold
        [~,c]=max(abs(mydir)); %Index of the maximum in absolute value
    else
        c=n+1;
    end
    if sign(mysign(c)) ~= sign(mydir(c))
        mydir=-mydir/abs(mydir(c));
    else
        mydir=mydir/abs(mydir(c));
    end
    dEst=mydir;
end

%Do one final fixed point correction after stepping to lambda=0
numfixedsteps=0;
error=inf;
while error>mythreshold && numfixedsteps<totsteps/10^3
        numfixedsteps=numfixedsteps+1;
        x=myEst(1:n);
        lambda=myEst(n+1);
        Rx=Rfun(F,x);
        Hess=3*Rx-lambda*Rz-(1-lambda)*Rzero;
        gradF=Rx*x-((lambda*Rz)+(1-lambda)*Rzero)*x;
        error=norm(gradF);
        myEst(1:n)=myEst(1:n)-Hess\gradF;
    end

myEstimate=myEst(1:n);
%Print lambda to ensure lambda=0, the number of steps used, and the norm of
%the gradient at the xestimate
lambda
numsteps
norm((Rfun(F,myEstimate)-Rzero)*myEstimate) %Gradient
end

%Below is some code to test the quadratic error criterion

% test=myEstimate;
% Rz=Rfun(F,z,5,20);
% %Create Rzero
% myZ = F.*(ones(5,1)*(y'));
% Rzero = (myZ*F')/20;
% Rx=Rfun(F,test,5,20);
% gradF=Rx*test-(Rzero)*test;
% norm(gradF)

%1/4<Rx*x,x>-1/2<R0x,x> + 1/4m
%(1/4)*dot(Rx*test,test)-(1/2)*dot(Rzero*test,test)+(1/(4*m))*sum(y.^2)


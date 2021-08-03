function [w,wmax,myC,condd] = solveForCredReal(Hess,Q,x,param,debug1)
%Solve [H,Qx]v=0 (v in null space)
%Choose c to be largest component in v, but marching in lambda if possible
%Hess is the Hessian Matrix, H(x)=3R(x)+lambda*Q-R0
%Q is an optional positve semi-definite matrix, nxn. Default to identity.
%x is is the signal of size n
%myV=abs(null([Hess,Q*x]));
n=param.n;
He=[Hess, Q*x];

if debug1.debugmode==1
    largeSv=svds(He,1);
    smallSv=svds(He,1,'smallest');
    condExt=cond(He);
end
[myQ,~,~]=qr(He','vector');
w=myQ(:,n+1);
[wmax,indexmax]=max(abs(w));
w=w/wmax;
myC=indexmax(1);
myCond=cond(Hess);
if myCond<param.cutOff
    myC=n+1;
end
if debug1.debugmode==1
    condd.myCond=myCond;
    condd.largeSv=largeSv;
    condd.smallSv=smallSv;
    condd.CondNumExt = condExt;
else
    condd=1;
end


% if cond(Hess)<param.cutOff
%     eta=-Hess\Q*x;
%     %Should satisfy the equation norm(myHess*eta+Q*x)=0, Hess*eta=-Qxiest
%     %-dot(Q*x,eta) this is the component on the 2n+1X2n+1 component of the hessian
%     w=[eta;1]; %This is the other null vector of the extended Hessian!
%     [wmax,indexmax]=max(abs(w));
%     w=(w/wmax);
%     myC=n+1;
% else %we march along a direction perpendicular to the lambda axis, according to the other vector in the null space
%     [Qred,Rred,~]=qr(Hess','vector');
%     zeta=Qred(1:n,n-1);
%     w=[zeta;0];
%     [wmax,indexmax]=max(abs(w));
%     w=(w/wmax);
%     size(w)
%     myC=indexmax(1);
end
 
 
 



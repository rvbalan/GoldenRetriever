function [xest,lambdaest,dxsign,dlambda_dt,trigger,debug1,debug2,param] = stepFullFun(F,Rzero,Q,param,xest,lambdaest,dxsign,dlambda_dt,y,debug1,debug2,essen)
%y here is a m-vector and is given by the CREATEy function.
%F is a matrix where the kth column represents the kth frame.
%Q is an optional positve semi-definite matrix, nxn
%deltat is a small positive real number which determines the step size
% The goal of this function is to follow the path of the derivatve being zero,
%taking a time step in (x,lambda), and then correcting the path
%cutOff will determine the cut off for numerical invertibility of a matrix
%The tolerance determines the convergence of the fixed point correction.
Myxold=xest;
Mylambdaold=lambdaest;


tolerance=param.tolerance;
cutOff=param.cutOff;
deltat=param.deltat;
n=param.n;
m=param.m;

if ~exist('Q','var')
    % Q parameter does not exist, so default it to the identity
    Q = eye(n);
end

%R0 = 1/m * sum(yk*fk*fk')

%Move (x,lambda) a time step.
myHess=Hess(F,xest,lambdaest,Q,Rzero,param);
if param.Certifier == 1
    [w,wmax,c,condd]=solveForCredRealcert(myHess,Q,xest,param);    %This solves [H,Qx]w=0 (w in null space)
else
    [w,wmax,c,condd]=solveForCredReal(myHess,Q,xest,param,debug1);    %This solves [H,Qx]w=0 (w in null space)
end

if debug1.debugmode==1
    debug1.CondNum=[debug1.CondNum, condd.myCond]; %Condition Number of Hessian
    debug1.extendedLargeSv=[debug1.extendedLargeSv, condd.largeSv]; % largest singular value
    debug1.extendedSmallSv=[debug1.extendedSmallSv, condd.smallSv]; % largest singular value
    debug1.CondNumExt=[debug1.CondNumExt, condd.CondNumExt]; 
    debug1.wnp1=[debug1.wnp1, w(n+1)];
    debug1.wc = [debug1.wc, w(c)];
end

if param.printIter==1
    myPercent=(lambdaest/eigs(Rzero,1))*100;
    if myPercent<0
        myPercent=0;
    end
    if myPercent>100
        myPercent=100;
    end
    Percent_Left=(100 - myPercent)
end

%fprintf(1,'\n Percent left %d \n',100 - myPercent)
%waitbar((100 - myPercent)/100,param.f,'Please wait...');

%if cond(myHess)<cutOff %test numerical invertibility with conditioning number
%myV=abs(null([myHess,Q*xest]));
%cutOff=-inf;
%if cond(myHess)<cutOff
if c==n+1 %Marched along lambda
    if param.printIter==1
        fprintf(1,'\n Marched along lambda direction \n')  
    end
    dxc_dt = sign(dlambda_dt)*sign(w(c)); %Unify this with dxsign by n+1 vector
else 
  %trigger.trigger=1; %Pass this information back to the runFunction so that it can do a different fixed point correction
   dxc_dt = sign(dxsign(c))*sign(w(c));
   if param.printIter==1
        fprintf(1,'\n Switched column number %d \n',c)
   end
end
    trigger.c=c; %Which component moving along
   %     if  backtrack.backtrackingmode==1
%         param.type{param.currectcounterstore+1}='Switched';
%         param.xstore{param.currectcounterstore+1}=xest;
%         param.lambdastore{param.currectcounterstore+1}=lambdaest;
%         param.currentcounterstore=param.currentcounterstore+1;
%     end
    
    %Chooses c to be largest component in v
    %stepsize=(param.arclength/norm(w));
    stepsize=(param.arclength/(abs(w(c))));
    %param.arclength
    stepsize=stepsize/(2)^(param.adjust);
    
    if lambdaest+w(n+1)*stepsize*dxc_dt<0
        stepsize=-lambdaest/(w(n+1)*dxc_dt);
        param.lambdaflag=1;
    end
    if debug1.debugmode==1
        debug1.dxsign=[debug1.dxsign,dxc_dt];
    end
    xest=xest+w(1:n)*stepsize*dxc_dt;
    lambdaest=lambdaest+w(n+1)*stepsize*dxc_dt;
    %dlambda_dt=sign(dxc_dt*sign(w(n+1)));
    dlambda_dt=sign(lambdaest-Mylambdaold);
    dxsign=sign(xest-Myxold);
    if  ~isempty(dxsign(dxsign==0))
        if param.printIter==1
            disp('Some coordinate does not change')
        end
        dxsign(dxsign==0)=1;
    end
%mydiffX=norm(Myxold-xest)
%mydiffL=norm(Mylambdaold-lambdaest)
%fprintf(debug1.fid,'%16.14f   %16.14f   %16.8f \n',mydiffX,mydiffL,wmax);
end





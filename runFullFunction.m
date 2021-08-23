function [soln,debug1,debug2,param,backtrack,debug3] = runFullFunction(y,F,param,essen,debug1,debug2,backtrack,debug3)
%y= vector of size m, create in the function createY()

%F= Matrix of the frame vectors as the columns, of dim n*m, created in
%MAIN()
%param = structure with the following componenets (all defined in MAIN):
%param.n= size n
%param.m= size m
%param.iterate = number of steps used in algorithm. Will divide lambdamax/param.iterate
%param.counter= Max number iterations allowed for fixed point correction
%param.cutoff=  Cutoff for condition number to switch paramaters (typically 10^4 or inf)
%param.tolerance = The error threshold in fixed point correction between successive iterations
%param.deltat = update size for xnew=xold + dx * deltat (same for lambda)
%param.Q = Q matrix used, usually the identity
%Also param stores debug mode cells which store what happened if something
%goes wrong.
%Write a script that calls this as a function to iterate over
%Insert debug paramater
% generate abs(F' * xest).^2 at different iterations, norm(error(y))
%To generate random Q, V=randn(n,2n). Q=(n/trace(V*V'))*V*V'
if debug1.debugmode==1
    debug1.eigHess1=[];
    debug1.eigHess2=[];
    debug1.eigHess3=[];
    debug1.eigHess4=[];
    debug1.eigHess5=[];
    debug1.wc=[];
    debug1.dxsign=[];
end

if ~exist('essen','var')
    % third parameter does not exist, so default it to something
    essen = 0;
end

Q=param.Q;
m=param.m;
n=param.n;

%Create Rzero
Z = F.*(ones(n,1)*(y'));
Rzero = (Z*F')/m;
if backtrack.backtrackingmode==1
    backtrack.Rzerostore{backtrack.loopIndex} = Rzero;
end


%Initialize lambda by lambda = eig max(Q^(-1/2)*R0*Q^(-1/2))
%C1=sqrtm(Q);
%C2=inv(C1);

%[emax,lambdamax]=eigs(C2 * Rzero * C2,1);
%emax=emax/norm(emax);

if param.eigStarter==1
    [emax,lambdamax]=eigs(Q\Rzero,1,'largestreal'); %Highest Eigenvalue
    %[lambdamax,indexmax]=max(diag(lambdamax));
    %emax=emax(:,indexmax);
else
    [emax,lambdamax]=eigs(Q\Rzero,param.eigStarter,'largestreal');
    %emax=emax(:,param.eigStarter);
    %lambdamax=lambdamax(param.eigStarter,param.eigStarter);
end

if backtrack.backtrackingmode==1
    backtrack.lm=lambdamax;
end



emax=emax/norm(emax);

%Q^(1/2)R0Q^(-1/2)*v1=lambda1*v1
%v1=Q^(1/2)*e1
%(lambda1*Q-R0)*e1=0
%[lambda1*I-Q^(-1)*R0]*e1=0

param.deltat=lambdamax/param.iterate;

%Initialize x
denom=0;
for k = 1:param.m
    denom = denom + abs(dot(emax,F(:,k)))^4;
end
denom = denom/param.m;
constantc=sqrt(param.deltat*dot(Q*emax,emax)/denom);
xest = constantc*emax;

% Initialize lambda
lambdaest = lambdamax-param.deltat;

%Fixed point correction for initialization
xest=fixedPointLambda(xest,lambdaest,param,F,Rzero);

%If the time step is 1
dxsign=sign(xest);
%Rename dlambda_dt
dlambda_dt=-1;


%Original Step Direction
Hessian=Hess(F,xest,lambdaest,Q,Rzero,param);
[~,wmaxzero,~] = solveForCredReal(Hessian,Q,xest,param,debug1);
param.wmaxzero=wmaxzero;

param.arclength=sqrt(param.deltat^2 + norm(xest)^2);
%param.arclength=param.deltat;
myTime=0; %Used to estimate time
myAlgCounter=0;
%for i=1:floor(4.6*(param.iterate-1)/5)
%param.f = waitbar(0,'1','Name','Percent Close to Lambda=0...');
param.firstFlag=1;
while lambdaest>0 && myAlgCounter<param.myCountThreshold && norm(xest)>(10^(-10)*sqrt(n))
%while lambdaest>0 && myAlgCounter<param.myCountThreshold

    myAlgCounter=myAlgCounter+1;
    if param.printIter==1
        fprintf('Your iteration number is %i', myAlgCounter)
    end
    i=myAlgCounter;
    %tic

    %If you want to run the certifier
    if param.Certifier == 0
        %First the step
        param.adjust = 0;
        [xest,lambdaest,dxsign,dlambda_dt,trigger,debug1,debug2,param]=stepFullFun(F,Rzero,Q,param,xest,lambdaest,dxsign,dlambda_dt,y,debug1,debug2);
        %Now Fixed Point Correction
        if param.Wolfe_Newton == 0
            [xest,lambdaest,myCount] = fixedPointAny(xest,lambdaest,param,F,Rzero,trigger.c);
        else
            [xest,lambdaest,myCount] = fixedPointNewton(xest,lambdaest,param,F,Rzero);
        end
    end
    
    
    
    
    if param.Certifier == 1 %If certifier is running, first compute things at xold above
        nx0 = norm(xest);
        %Can compute these once in the beginning
        Fcol=vecnorm(F).^2; %Norm of columns of F
        B = norm(F,2)^2; %Upper bound on B, the frame bound
        beta=max(Fcol)*B; %Beta in the notes %Upper bound on B=max_k||f_k||
        Hess0 = Hess(F,xest,lambdaest,Q,Rzero,param);
        Hessext = [Hess0,Q*xest];
        smin =svds(Hessext,1,'smallest');
        S0=svds(Hessext,1); %||H_{ext,0}||
        failedflag=0;
        oldarclength = param.arclength;
        oldx=xest;
        oldlambda = lambdaest;
        for stepCert=1:10
            param.adjust = 0;
            %param.adjust = stepCert; old case, now we divide the arclength
            %into 2
            if stepCert == 10 %If reach 10 steps, break
                failedflag=1;
                break
            end
            %Save the old variables to certify path later
            %First the step
            
            [xest,lambdaest,dxsign,dlambda_dt,trigger,debug1,debug2,param]=...
                stepFullFun(F,Rzero,Q,param,oldx,oldlambda,dxsign,dlambda_dt,y,debug1,debug2);
            
            %Compute the reduced Hessian using the column at the step
            
            Hessred=Hessext;
            Hessred(:,trigger.c) = []; %Gives H_{red:c,0}
            sminred = svds(Hessred,1,'smallest');
            
            
            
            %Now Fixed Point Correction
            if param.Wolfe_Newton == 0
                [xest,lambdaest,myCount] = fixedPointAny(xest,lambdaest,param,F,Rzero,trigger.c);
            else
                [xest,lambdaest,myCount] = fixedPointNewton(xest,lambdaest,param,F,Rzero);
            end
            [certificate,newstep,param] = ...
                CertifierFullA(xest,lambdaest,oldx,oldlambda,F,smin,S0,Q,trigger.c,Hessext,param);
            
            if certificate == 1
                %stepCert
                break
            else
                param.arclength = (param.arclength/(4)); %Half the step size
                %dnd=param.arclength
            end
        end
        param.adjust = 1;
        param.arclength=oldarclength;
        if failedflag == 1
            failedflag
        end
    end %END OF CERTIFIER
    
    
    %Save Debug and Backtrack Information
    if debug1.debugmode1==1
        myHess=Hess(F,xest,lambdaest,Q,Rzero,param);
        [myNorm,myIX,mygradJ,myJX]=calcDebuger1(xest,y,F,param,Rzero,lambdaest);
        debug1.normX=[debug1.normX myNorm]; %Tracks the norm of the estimate x^t
        debug1.Ix=[debug1.Ix myIX]; %Tracks value of I(x)=1/4m * sum(|yk-|<x,fk>|^2|^2)
        debug1.gradJ=[debug1.gradJ mygradJ]; %Tracks the norm of the gradient
        debug1.Jx=[debug1.Jx myJX];
        debug1.eig=[debug1.eig eigs(myHess,1,'smallestreal')];%Tracks the min eig of Hessian
        %Comment out this section
%         A=sort(eig(myHess),'descend');
%         debug1.eigHess1=[debug1.eigHess1, A(1)];
%         debug1.eigHess2=[debug1.eigHess2, A(2)];
%         debug1.eigHess3=[debug1.eigHess3, A(3)];
%         debug1.eigHess4=[debug1.eigHess4, A(4)];
%         debug1.eigHess5=[debug1.eigHess5, A(5)];
%         
        
        
        debug1.fixedPointCounter=[debug1.fixedPointCounter myCount];
        debug1.lambdaValue=[debug1.lambdaValue lambdaest];
        debug1.Cstore=[debug1.Cstore trigger.c];
        
        if myCount>=param.counter && backtrack.backtrackingmode==1
            backtrack.type{backtrack.currentcounterstore+1}='Fixed Point not Converge';
            backtrack.xstore{backtrack.currentcounterstore+1}=xest;
            backtrack.lambdastore{backtrack.currentcounterstore+1}=lambdaest;
            backtrack.iterationNum{backtrack.currentcounterstore+1} = i;
            backtrack.tracker{backtrack.currentcounterstore+1} = backtrack.loopIndex;
            backtrack.currentcounterstore=backtrack.currentcounterstore+1;
        end
        if eigs(myHess,1,'smallestreal')<=0 && backtrack.backtrackingmode==1
            backtrack.type{backtrack.currentcounterstore+1}='Minimum Eig non-positive';
            backtrack.xstore{backtrack.currentcounterstore+1}=xest;
            backtrack.lambdastore{backtrack.currentcounterstore+1}=lambdaest;
            backtrack.iterationNum{backtrack.currentcounterstore+1} = i;
            backtrack.tracker{backtrack.currentcounterstore+1} = backtrack.loopIndex;
            backtrack.currentcounterstore=backtrack.currentcounterstore+1;
        end
        
    end
    if debug2.debugmode2==1
        [myMinx]=calcDebuger2(xest,param,lambdaest,essen);
        debug2.minx=[debug2.minx myMinx];
        
        %Method 1
        %t=(3*norm(xest)^2 + lambdaest - lambdamax)/(3*norm(essen.x)^2 - lambdamax);%Tracking the M(t) matrix path... Approach 1
        %if t<0
        %    t=0;
        %end
        %if t>1
        %    t=1;
        %end
        % End Method 1
        %Method 2
        %         if i==1
        %             t=0;
        %         else
        %             t=tIterate(t,norm(xest)^(2),lambdaest,norm(essen.x)^(2),Rzero,essen.x);
        %         end
        %         %End Method 2
        %          myz=essen.x;
        %          M=(1-t)*Rzero + 2*t*essen.x*essen.x';
        %          mu=eigs(M,1);
        %          eigmu=null(M-mu*eye(size(M,1)));
        %          xt=(sqrt(t*norm(myz)^2 - lambdaest+mu)/sqrt(3))*eigmu;
        %          closert=min(norm(xest-xt)^2,norm(xest+xt)^2);
        %          debug2.diffxt=[debug2.diffxt closert];
        %          debug2.mytValue=[debug2.mytValue t];
    end
    if debug3.debugmode3==1
        debug3.mySavedX{i}=xest;
    end
    
    %     myTime=myTime+toc;
    %     TimeLeft=((param.myCountThreshold/i)-1)*myTime;
    %     HoursLeft=floor(TimeLeft/3600);
    %     MinutesLeft=floor(60*rem(TimeLeft/3600,1));
    %     SecondsLeft= TimeLeft-(3600*HoursLeft+60*MinutesLeft);
    %fprintf('Your Estimated Time Left is %d Hours and %d Minutes and %f Seconds !\n',HoursLeft,MinutesLeft,round(SecondsLeft*100)/100)
end
%Get eigenvalue information
myHess=Hess(F,xest,lambdaest,Q,Rzero,param);
param.Heig=eigs(myHess,1,'smallestreal');
%Solution
soln.x=xest;
soln.lambda=lambdaest;
soln.deltat=param.deltat;
end

function [err,xestimate,backtrack,debug3,debug1] = Recoverer(y,F,x,snr,Q,loopIndex,currentCounter)
%This Function solves for a given y input as well as the frame F. The x and
%snr inputs are optional for additional debugging

%loopIndex is if you are interating over the Recoverer, you can save your
%backtracking information (frame/Rzero/Q) with the index you use. It
%differs from currentCounter because it only saves things that are
%initialized once in the beginning of the algorithm, and is not updated
%until Recoverer is called again.

%currentCounter is used in backtracking to continue where you left off, if
%running under multiple loops. This iterates every time something comes up
%for reproducability (like the current x, current lambda,iteration number,
%etc.) It does not include the loopIndex (frame/Rzero/Q), because those
%only change upon different input, not different loops within the code.
%Upon calling the Recoverer multiple times, set currentCounter=0, then set
%currentCounter=backtrack.currentCounter and pass that in for successive
%loops.


%Most important Parameters for changing
param.iterate=1e5; %Used to determine arclength, divide lambdamax/param.iterate
param.eigStarter=1; %Which Eigenvalue you want to start from (default is 1)
param.printIter=0; %Print output during iterations? 0 for no, 1 for yes
param.Wolfe_Newton = 1; %0 for Frank Wolfe ... 1 for Newton, best to use Newton
param.lambdaflag=0; %This identifies whether lambda ever achieves lambda=0
param.WF=0; %Do you want to run Wirtinger Flow after you reach lambda=0
debug1.debugmode=0; %Do you want to run with debug and backtrack (1 for active, 0 for inactive)
%To plot the difference, need to set debug1.debugmode to 1


%Often don't use loopIndex or currentCounter unless they are in one giant
%loop, so safe to ignore this
if ~exist('loopIndex','var')
  loopIndex=1;
end
if ~exist('currentCounter','var')
  currentCounter=0;
end
if ~exist('Q','var')
  Q=eye(size(F,1));
end


param.Certifier=0; %Do not use yet, still needs work

%Initialize other debug modes, later version will change how this is done
if debug1.debugmode==1
    debug1.debugmode1=1;
    debug2.debugmode2=1;
    debug3.debugmode3=1;
    backtrack.backtrackingmode=1;
else
    debug1.debugmode1=0;
    debug2.debugmode2=0;
    debug3.debugmode3=0;
    backtrack.backtrackingmode=0;
end

%If no initial x or snr is given, don't crash the program.

if ~exist('x','var')
  x=0;
  myTriggerX=0;
  debug2.debugmode2=0;
end
if ~exist('snr','var')
  snr=inf;
  myTriggerSNR=0;
end

%Set the less important paramters of the algorithm

param.n=size(F,1); % Set size of the vectors,
param.m=size(F,2);
param.counter=50; %Max number needed for fixed point correction
param.cutOff=50; %Typically 50, if you want to switch to a different direction (or 10^4 for being more conservative)2
param.tolerance=10^-7; %The error threshold in fixed point correction between successive iterations
param.Q=Q; %The Q value used, usually the identity matrix
%param.myCountThreshold=param.iterate*4; % used to bound the total number of iterations before stopping
param.myCountThreshold=param.iterate*4;
param.rho = 1/2;
%essen stores the true solution and the snr, something a blind solver
%wouldn't have access to.
essen.x=x;
essen.snr=snr;

%Backtrack stores all the data associated with a certain point in the
%iteration assuming certain criteria are triggered.

%Following functions initialize the backtracking and the debugging
%information for the algorithm

if backtrack.backtrackingmode==1
backtrack=initializeBacktrack(backtrack,currentCounter,loopIndex,F,param,y,x);
end
if debug1.debugmode==1
[debug1,debug2,debug3] = initializeDebug(debug1.debugmode1,debug2.debugmode2,debug3.debugmode3,backtrack,param,F,y,x); 
end


[soln,debug1,debug2,param,backtrack,debug3]=runFullFunction(y,F,param,essen,debug1,debug2,backtrack,debug3);
xestimate=soln.x;

%If you want to do Wirtinger Flow from the solution of the golden retriever
if param.WF==1
    if param.lambdaflag==1   
        normest = sqrt(sum(y)/numel(y));
        T = 25000;                           % Max number of iterations
        tau0 = 330;                         % Time constant for step size
        mu = @(t) min(1-exp(-t/tau0), 0.002); % Schedule for step size
        for t = 1:T
            yz = F'*xestimate;
            grad  = (1/param.m)* F*( ( abs(yz).^2-y ) .* yz ); % Wirtinger gradient
            xestimate = xestimate - mu(t)/normest^2 * grad;             % Gradient update
        end
    end
end
err=min(norm(xestimate-x),norm(xestimate+x));

%fclose(debug1.fid);
%Plots of the tracked information in the debug modes
%Uncomment
if debug1.debugmode==1
    Plotter(debug1,debug2,debug3) %Plot the solution
end
end
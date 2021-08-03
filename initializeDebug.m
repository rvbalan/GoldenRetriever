function [debug1,debug2,debug3] = initializeDebug(debugmode1,debugmode2,debugmode3,backtrack,param,F,y,x)
%Initializes all the debug information

%Debug 1 first (one that can always be tracked)
%debug1 can always be activated for tracking information
debug1.debugmode1=debugmode1;
debug1.debugmode=1;
debug1.eig=[];%Tracks the min eig of Hessian
debug1.gradJ=[]; %Tracks the norm of the gradient
debug1.normX=[]; %Tracks the norm of the estimate x^t
debug1.Ix=[]; %Tracks value of I(x)=1/4m * sum(|yk-|<x,fk>|^2|^2)
debug1.Jx=[]; %Tracks the J criterion, I(x) + lambda/2<Qx,x>
debug1.fixedPointCounter=[]; %Tracks the number of fixed point corrections used
debug1.lambdastore=[]; %Keeps track of lambda over time
debug1.Cstore=[]; %Keeps track of which column switched over time
debug1.extendedLargeSv=[]; % largest singular value
debug1.extendedSmallSv=[]; % largest singular value
debug1.CondNum=[]; %Condition Number of Hessian
debug1.CondNumExt=[];
debug1.lambdaValue=[];
debug1.wnp1=[]; % Lambda coordinate w(n+1) in null vector

%Debug 2 (if true solution is known) 
%debug2 "cheats" and gets additional information about the true solution
%Only can be used if true solution is known
debug2.debugmode2=debugmode2;
debug2.minx=[]; %Tracks the difference between estimate and true x

%Debug 3 (for memory intensive investigation)
%ebug3 is saving the entire path in the x-space. It can always be used, but takes up more space.
debug3.debugmode3=debugmode3;
debug3.mySavedX=cell(1,1);

%Backtrack
%backtrack.currentcounterstore=currentCounter;
% backtrack.xstore=cell(1,1);
% backtrack.framestore=cell(1,1);
% backtrack.lambdastore{backtrack.currentcounterstore+1}=cell(1,1);
% backtrack.Qstore{backtrack.loopIndex}=cell(1,1);
% backtrack.Rzerostore{backtrack.loopIndex} = cell(1,1);
% backtrack.ystore{backtrack.loopIndex} = cell(1,1);
% backtrack.type{backtrack.currentcounterstore+1} = cell(1,1);
% backtrack.iterationNum{backtrack.currentcounterstore+1} = cell(1,1);
% backtrack.tracker{backtrack.currentcounterstore+1} = cell(1,1); %The tracker just saves which loopIndex is being used for each of the outputs.
% if backtrack.backtrackingmode==1
%     backtrack.framestore{backtrack.loopIndex}=F;
%     backtrack.Qstore{backtrack.loopIndex}=param.Q;
%     backtrack.ystore{backtrack.loopIndex}=y;
%     backtrack.zstore{backtrack.loopIndex}=x; %Clean Signal x, we save it under z to avoid confusion with any future x
% end



end


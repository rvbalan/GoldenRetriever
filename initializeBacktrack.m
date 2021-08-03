function [backtrack] = initializeBacktrack(backtrack,currentCounter,loopIndex,F,param,y,x)
%Initializes all the backtrack information
backtrack.currentcounterstore=currentCounter;
backtrack.loopIndex=loopIndex;
backtrack.xstore=cell(1,1);
backtrack.framestore=cell(1,1);
backtrack.lambdastore{backtrack.currentcounterstore+1}=cell(1,1);
backtrack.Qstore{backtrack.loopIndex}=cell(1,1);
backtrack.Rzerostore{backtrack.loopIndex} = cell(1,1);
backtrack.ystore{backtrack.loopIndex} = cell(1,1);
backtrack.type{backtrack.currentcounterstore+1} = cell(1,1);
backtrack.iterationNum{backtrack.currentcounterstore+1} = cell(1,1);
backtrack.tracker{backtrack.currentcounterstore+1} = cell(1,1); %The tracker just saves which loopIndex is being used for each of the outputs.
if backtrack.backtrackingmode==1
    backtrack.framestore{backtrack.loopIndex}=F;
    backtrack.Qstore{backtrack.loopIndex}=param.Q;
    backtrack.ystore{backtrack.loopIndex}=y;
    backtrack.zstore{backtrack.loopIndex}=x; %Clean Signal x, we save it under z to avoid confusion with any future x
end
end


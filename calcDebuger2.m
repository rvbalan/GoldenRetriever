function [myMinx]=calcDebuger2(xest,param,lambdaest,essen)

myMinx=min([norm(xest-essen.x)^2,norm(essen.x+xest)^2]);

end



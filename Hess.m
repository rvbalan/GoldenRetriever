function [myHess] = Hess(F,xest,lambda,Q,Rzero,param)
%This will compute the Hessian at a given (x,lambda) value
%F is a matrix of the frame vectors as the columns, of dim n*m, created in
%MAIN()
%xest is the vector input (size n) which you want to compute the hessian
%at.
%lambda is a single number which is input into the Hessian
%Q is a positive definite symmetric matrix, which is usually default to the
%identity.
%Rzero is a matrix (nxn) which is defined in runFunction()
%param is a structure which contains the sizes, m and n, among other
%parameters
myHess=3*Rfun(F,xest,param) + lambda*Q - Rzero;
end
 



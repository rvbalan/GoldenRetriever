function [Rx] = Rfun(F,xest,param)
%Here we create R(xest)= 1/m * sum(|<xest,fk>|^2 * fk* fk')
%F is the matrix of the frame vectors as the columns, of dim n*m, created in
%MAIN()
%xest is the vector input in function above, of size n
%param is a structure parameter, which contains the vector size n, and the
%frame size m, among other things.
if ~exist('param','var')
  [n,m]=size(F);
else
  m=param.m;
  n=param.n;
end

cx = (abs(F'*xest)).^2; % Coefficients squared, i.e. an m vector with <xest,f_k>^2
%Rx = (F*diag(cx)*F')/m; %=1/m \sum <xest,f_k>^2 f_k f_k'

Z = F.*(ones(n,1)*(cx'));
Rx = (Z*F')/m;
% Z=zeros(m,n);
% for i=1:m
%     Z(i,:)=(cx(i)*F(:,i))';
% end
% Rx = (F*Z)/m; %=1/m \sum <xest,f_k>^2 f_k f_k'
% end




%function [Rx] = Rfun(F,xest,param)
%%Here we create R(xest)= 1/m * sum(|<xest,fk>|^2 * fk* fk')
%m=param.m;
%n=param.n;
%Rx=zeros(n); %Initialize R
%for k =1:m
 %   Rx = Rx + ((F(1:n,k)' * xest)^2)*F(1:n,k)*(F(1:n,k))';
%end
%Rx = Rx/m;
%end



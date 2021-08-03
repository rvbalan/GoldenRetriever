function [] = plotEnd(truez,F,myy,lm,Q)
%rho1color=[.949,.027,.074];
%GRcolor = [.984,.603,.054];
GRcolor=[.949,.027,.074];
rho1color = [.984,.813,.039];
myview=[.8,-.1,.1];
figure(1);
view(myview)
n=size(truez,1);
m=size(F,2);
%First plot the plane and the line
[y, z] = meshgrid(-1.5:0.1:1.5); % Generate x and y data
x = zeros(size(y, 1)); % Generate z data
%surf(x, y, z,'Edgecolor','none') % Plot the surface
surf(x, y, z,'Edgecolor',[1,1,1]) % Plot the surface
colormap(winter)
alpha(.5)
hold on
xt = @(t) t;
yt = @(t) 0;
zt = @(t) 0;
hold on
essen.x=truez;
essen.snr=inf;
if ~exist('myy','var')
  [myy,~] = createYdb(essen,F);
end
if ~exist('Q','var')
  Q=eye(length(truez));
end
[err,~,backtrack,debug3,debug1] = Recoverer(myy,F,truez,inf,Q);
myerr=err

if ~exist('lm','var')
  lm=backtrack.lm;
end
%plot the line from 0 to lm
fplot3(xt,yt,zt,[0,lm],"k","linewidth",3)
hold on

% lparam= lm-8.5988*(y.^2+z.^2);
% surf(lparam,y,z,'facealpha',.1)
% 
% end

AllX=zeros(length(debug3.mySavedX),n);
for tempindex=1:length(debug3.mySavedX)
    AllX(tempindex,:) = debug3.mySavedX{tempindex};
end
[USVD,S,VSVD]=svd(AllX,'econ');
myMatrix=AllX*VSVD;
mye1=myMatrix(1,:)';
mye2=myMatrix(2,:)';

%Now we plot the curve in gold
previous=[lm,0,0];
length(debug3.mySavedX)
interval=floor(length(debug3.mySavedX)/90);
interval
if interval==0
    interval=1;
end

for i=1:interval:length(debug3.mySavedX)
    backSol=debug3.mySavedX{i};
    Bs1=dot(backSol,mye1);
    Bs2=dot(backSol,mye2);
    lambda=debug1.lambdaValue(i);
    zerot=@(t) 0;
    xt = @(t) lambda-t*lambda+t*previous(1);
    %yt = @(t) backSol(1)-t*backSol(1)+t*previous(2);
    %zt = @(t) backSol(2)-t*backSol(2)+t*previous(3);
    yt = @(t) Bs1-t*Bs1+t*previous(2);
    zt = @(t) Bs2-t*Bs2+t*previous(3);
    %plot the line from between consecutive points.
    fplot3(xt,yt,zt,[0,1],'Color',GRcolor,"linewidth",1)
    plot3(lambda,Bs1,Bs2,'Color',GRcolor,'Marker','.','MarkerSize',10)
    fplot3(zerot,yt,zt,[0,1],'Color','g',"linewidth",1)
    plot3(0,Bs1,Bs2,'Color','g','Marker','.','MarkerSize',10)
    hold on
    previous=[lambda,Bs1,Bs2];
end
view(myview)
end


function [] = Plotter(debug1,debug2,debug3)

%vec=eigs(debug3.Rzero);

if debug1.debugmode1==1
    figure(2)
    subplot(4,3,1)
    plot(debug1.eig)
    title('Minimum eig of Hessian')
    subplot(4,3,2)
    plot(debug1.Ix)
    title('Value of I(x) criterion')
    subplot(4,3,3)
    plot(debug1.gradJ)
    title('Value of norm(gradJ)')
    subplot(4,3,4)
    plot(debug1.normX)
    title('Value of norm(x)^2')
    subplot(4,3,5)
    plot(debug1.Jx)
    title('Value of J(x) criterion')
    subplot(4,3,6)
    plot(debug1.fixedPointCounter)
    title('Number of Fixed Points used')
    subplot(4,3,7)
    plot(debug1.lambdaValue)
    title('Lambda Value')
%     hold on
%     plot(0,vec(2),'r*')
%     hold on
%     plot(0,vec(3),'r*')
%     hold on
%     plot(0,vec(4),'r*')
%     hold on
%     plot(0,vec(5),'r*')
    subplot(4,3,8)
    plot(debug1.Cstore)
    title('Column which switched')
    subplot(4,3,9)
    plot(debug1.wnp1)
    title('w(n+1)')
    subplot(4,3,10)
    plot(debug1.CondNum)
    title('Cond number of Hessian')
    subplot(4,3,11)
    plot(debug1.CondNumExt)
    title('Cond number of Extended Hessian')
end
if debug2.debugmode2==1
    subplot(4,3,12)
    plot(debug2.minx)
    title('Value of difference between (estimate and true x)^2')
end


figure(3)
plot(debug1.eigHess1,'DisplayName','Eig1')
hold on
plot(debug1.eigHess2,'DisplayName','Eig2')
plot(debug1.eigHess3,'DisplayName','Eig3')
plot(debug1.eigHess4,'DisplayName','Eig4')
plot(debug1.eigHess5,'DisplayName','Eig5')
title('Eigenvalues of Hessian')
hold off
lgd = legend;

figure(4)
plot(debug1.wc)
title('wc')

figure(5)
plot(debug1.dxsign)
title('dxsign')


end


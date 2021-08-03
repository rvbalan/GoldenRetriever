function [err,x,F,y] = GR_fun(n,m,x,F,snr,Q,y)
%This generates a random instance of x and F with given snr and runs it

if ~exist('Q','var')
  Q=eye(n);
end
if ~exist('snr','var')
  essen.snr=inf;
else
    essen.snr=snr;
end
if ~exist('F','var')
  F=randn(n,m);
end
if ~exist('x','var')
  x=1/sqrt(n) * randn(n,1);
end
if x=='no'
    x=1/sqrt(n) * randn(n,1);
end
if ~exist('y','var')
    essen.x=x;
    [y,mysigma] = createYdb(essen,F);
else
    essen.x=x;
end




err = Recoverer(y,F,essen.x,essen.snr,Q);
%[~,errorfactor] = computePhase(xestimate,essen.x);
%errorfactor=min(norm(xestimate+essen.x),norm(xestimate-essen.x));
%errorfactor=errorfactor^2;

end


function [y,mysigma] = createYdb(essen,F)
%Input format is (essen,F,parm)
%essen is a structure which contains the snr, as well as
%x which is the desired signal 
%F is the nxm matrix where the columns are fk, the frame vectors
%param is a structure which contains the sizes, m and n, among other
%parameters
snr=essen.snr;
m=size(F,2);
x=essen.x;
%rng(essen.seed);
noise=randn(m,1);
 
magSquare = abs(F' * x).^2;
 
%sigma is the level of noise you want, computed from snr
if (snr == inf)
    sigma = 0;
else
    % snrdb = 20*log_10(norm(magSquare) / (sigma*norm(noise)))
            sigma = norm(magSquare) / (10^(snr/20)*norm(noise));
end
mysigma=sigma;
y = magSquare + sigma*noise;
end
 
 
 
 




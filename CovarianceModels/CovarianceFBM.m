function C = CovarianceFBM(N,D,alpha,sigma,R,dt)

if R == 0 % covariance matrix without localization noise
    C = zeros(N,N);
    model = 2*Dsim;
    for n = 1:N
        for m = 1:N
            C(n,m) = model/2*dt^alpha*(abs(n-m+1)^alpha - 2*abs(n-m)^alpha + abs(n-m-1)^alpha);
        end
    end
    
else % covariance matrix with localization noise
    x = 2:N-1;
    model = D*dt^alpha/((alpha+2)*(alpha+1));
    A = @(x,alpha) (x+1).^(alpha+2) + (x-1).^(alpha+2) - 2*x.^(alpha+2);
    C1 = 2*model*(A(1,alpha) - 2) + 2*sigma^2;
    C2 = model*(A(2,alpha)-2*A(1,alpha) + 2) - sigma^2;
    C3 = model*(A(1+x,alpha) - 2*A(x,alpha) + A(x-1,alpha));
    C = toeplitz([C1 C2 C3]);
end

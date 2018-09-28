function Ahat = LRcompress(A,tol,tolabs)   




[m,n] = size(A);


[U,S,V] = svd(A,0);
s = diag(S);
if nargin == 3
if(s(1)<tolabs)
    r=0;
    return
end
end

r = min(m,n);
for i=1:min(m,n)
    if (s(i)/s(1)<=tol)
        r=i;
        if(s(i)<s(1)*tol/10)
            r = i -1;
        end
        break
    end
end

Ahat = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';   


end
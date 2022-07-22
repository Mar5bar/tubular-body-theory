function [z2p] = Z2p(as,alphas,us)

z2 = Z2(as,alphas);

z2p = zeros(size(as));

for ind = 1 : numel(as)
    a = as(ind);
    alpha = alphas(ind);
    u = us(ind);
        
    k1 = sqrt(1 - u);
    k2 = sqrt(1 + u);
    
    if abs(alpha - 1) < 1e-6
        if abs(u) < 1e-6
            zt2=sqrt(2);
        elseif abs(u - 1) < 1e-6 || abs(u + 1) < 1e-6
            zt2=2;
        else
            zt2=sqrt(2)*(k2-k1)./u;
        end
    else
        
        l=log(-(((1+u-alpha^2+sqrt((1+u).*(1+u-2*u*alpha^2+(u-1).*alpha^4))).*(u-1+alpha^2-sqrt((u-1).*(alpha^4-1+u.*(alpha^2-1)^2))))./((alpha^2).*(2-2*alpha^2 + (u.^2).*(alpha^2-2)))));
        
        zt2=l./sqrt(1-alpha^2);
    end
    
    z2p(ind) =real(z2(ind)-2*pi*zt2/a);
    
end
end


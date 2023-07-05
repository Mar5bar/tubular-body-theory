function [z1p] = Z1p(as,alphas,us)

z1 = Z1(as,alphas);

z1p = zeros(size(as));

for ind = 1 : numel(as)
    a = as(ind);
    alpha = alphas(ind);
    u = us(ind);
        
    k1 = sqrt(1 - u);
    k2 = sqrt(1 + u);
    
    if abs(alpha - 1) < 1e-6

        if abs(u) < 1e-6
            zt1=7/(3*sqrt(2));
        elseif abs(u - 1) < 1e-6 || abs(u + 1) < 1e-6
            zt1=8/3;
        else
            zt1= (8*(k2-k1)+u.*(-4*(k1+k2)+u.*(19*(k1-k2)+3*u.*(u.*(k2-k1)+4*(k1+k2))))) ./ (3*(u.^3).*sqrt(2-2*u.^2));
        end

    else
        
        l=log(-(((1+u-alpha^2+sqrt((1+u).*(1+u-2*u*alpha^2+(u-1).*alpha^4))).*(u-1+alpha^2-sqrt((u-1).*(alpha^4-1+u.*(alpha^2-1)^2))))./((alpha^2).*(2-2*alpha^2 + (u.^2).*(alpha^2-2)))));
        
        if abs(u - 1) < 1e-6 || abs(u + 1) < 1e-6
            zt1=2/(-1+alpha^2)-(2*(-2+alpha^2)*log((1+sqrt(1-alpha^2))/alpha))/(1-alpha^2)^(3/2);
        else
            zt11=-((l.*(-2+alpha^2))/(1-alpha^2)^(3/2));
            zt12=-(2*k1*(-(k1^2)*(k2^4)+(1+2*u-u^3)*alpha^2)*sqrt(1-alpha^4-u*((-1+alpha^2)^2))+2*k2*(-(k1^4)*(k2^2)+(1-2*u+u^3)*alpha^2)*sqrt(1-alpha^4+u*((-1+alpha^2)^2)))/((-1+alpha^2)*(2-2*alpha^2+(u^2)*(-2+alpha^2))*sqrt((1-alpha^2)*((u^4)*((-1+alpha^2)^2)+((1+alpha^2)^2)-2*(u^2)*(1+alpha^4))));
            zt1=zt11+zt12;
        end

    end
    
    z1p(ind) = real(z1(ind)-2*pi*zt1/a);
end

end


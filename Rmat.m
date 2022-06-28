function [R, fs, fsTotal, S, Phi] = Rmat(SO0,SNO,N,numS,numPhi,epps,rho,r1,r2,r3,erho1,erho2,erho3)
%RMATV3 Summary of this function goes here
%   Detailed explanation goes here

R=zeros(6,6,N+1);

S1=@(s,phi)r1(s)+epps*rho(s).*erho1(s,phi);
S2=@(s,phi)r2(s)+epps*rho(s).*erho2(s,phi);
S3=@(s,phi)r3(s)+epps*rho(s).*erho3(s,phi);

ds = 2/(numS); %distance between consecutive s points.
s = linspace(-1+(ds/2),1-(ds/2),numS); %s points to be used for colocation method
dphi = 2*pi/(numPhi);
phi=linspace(-pi+dphi/2,pi-dphi/2,numPhi);
[S,Phi]=meshgrid(s,phi);

x = S1(S',Phi'); % The ordering in the matrices is the opposite to meshgrid.
y = S2(S',Phi');
z = S3(S',Phi');
X = [x(:)';y(:)';z(:)'];
XRep = repmat(X,1,1,6);

v1=@(s,phi) [1,0,0];
v2=@(s,phi) [0,1,0];
v3=@(s,phi) [0,0,1];
v4=@(s,phi) cross([1,0,0],[S1(s,phi),S2(s,phi),S3(s,phi)]);
v5=@(s,phi) cross([0,1,0],[S1(s,phi),S2(s,phi),S3(s,phi)]);
v6=@(s,phi) cross([0,0,1],[S1(s,phi),S2(s,phi),S3(s,phi)]);


V0 = zeros(3*(numPhi)*numS,6);
for i=1:numS
    for j=1:(numPhi)
        V0((1:3)+3*(i-1)+3*numS*(j-1),1)=8*pi*v1(s(i),phi(j));
        V0((1:3)+3*(i-1)+3*numS*(j-1),2)=8*pi*v2(s(i),phi(j));
        V0((1:3)+3*(i-1)+3*numS*(j-1),3)=8*pi*v3(s(i),phi(j));
        V0((1:3)+3*(i-1)+3*numS*(j-1),4)=8*pi*v4(s(i),phi(j));
        V0((1:3)+3*(i-1)+3*numS*(j-1),5)=8*pi*v5(s(i),phi(j));
        V0((1:3)+3*(i-1)+3*numS*(j-1),6)=8*pi*v6(s(i),phi(j));
    end
end

f0=SO0*V0;

fs = cell(N+1,1);
fs{1} = f0;
fsTotal = cell(N+1,1);
fsTotal{1} = f0;

% Compute R using f_0.
R(1:3,:,1) = squeeze(sum(reshape(f0,3,[],6),2)) * dphi * ds; % Total force.
R(4:6,:,1) = squeeze(sum(cross(XRep,reshape(f0,3,[],6)),2)) * dphi * ds; % Total torque.

% Iterate to construct R using f_n.
fo=f0;
sgn = -1;
for n = 1 : N
    fi = SNO*fo;
    fs{n+1} = fi;
    fsTotal{n+1} = fsTotal{n} + sgn * fi;
    
    R(1:3,:,n+1) = R(1:3,:,n) + sgn * squeeze(sum(reshape(fi,3,[],6),2)) * dphi * ds; % Total force.
    R(4:6,:,n+1) = R(4:6,:,n) + sgn * squeeze(sum(cross(XRep,reshape(fi,3,[],6)),2)) * dphi * ds; % Total torque.
    
    fo=fi;
    sgn = -sgn;
end

end


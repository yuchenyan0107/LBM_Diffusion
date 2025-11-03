clc
close all
clear all

syms rho ux uy u vx vy v lam
phi=sym('phi','positive');
u=[ux;uy];
v=[vx;vy];

fm=rho/(2*pi*phi/3)*exp(-(3*sum((v-u).^2))/(2*phi));
m0  =simple(int(int( fm ,vx,-inf,inf),vy,-inf,inf))

m1x =simple(int(int( vx*fm ,vx,-inf,inf),vy,-inf,inf))
m1y =simple(int(int( vy*fm ,vx,-inf,inf),vy,-inf,inf))

m2xx=simple(int(int( vx*vx*fm ,vx,-inf,inf),vy,-inf,inf))
m2yy=simple(int(int( vy*vy*fm ,vx,-inf,inf),vy,-inf,inf))
m2xy=simple(int(int( vx*vy*fm ,vx,-inf,inf),vy,-inf,inf))

m3x =simple(int(int( (vy^2)*vx*fm ,vx,-inf,inf),vy,-inf,inf))
m3y =simple(int(int( (vx^2)*vy*fm ,vx,-inf,inf),vy,-inf,inf))
m4  =simple(int(int( (vx^2)*(vy^2)*fm ,vx,-inf,inf),vy,-inf,inf))

% Original
%=========================
% mx3 = 1/3*ux*rho*(3*uy^2+phi)
% m3y = 1/3*uy*rho*(3*ux^2+phi)

% Simplified
%=========================
% cutting 3rd order terms in ux and uy and phi = 1
m3x = 1/3*ux*rho;
m3y = 1/3*uy*rho;

% Original
%=========================
% m4  = phi*(ux^2*uy^2*rho+1/3*ux^2*rho+1/3*uy^2*rho+1/9*rho*phi)

% Simplified
%=========================
% cutting 4th order terms in ux and uy and common phi = 1
m4  = 1/3*(ux^2+uy^2)*rho+1/9*rho*phi;

m=[m0;m1x;m1y;m2xx;m2yy;m2xy;m3x;m3y;m4];

V=sym([0     0
     1     0
     0     1
    -1     0
     0    -1
     1     1
    -1     1
    -1    -1
     1    -1]);

M = [ones(1,9);...
    V(:,1)';...
    V(:,2)';...
    (V(:,1).*V(:,1))';...
    (V(:,2).*V(:,2))';...
    (V(:,1).*V(:,2))';...
    (V(:,2).^2)'.*V(:,1)';...
    (V(:,1).^2)'.*V(:,2)';...
    (V(:,1).^2)'.*(V(:,2).^2)'];

fmd = simple(expand(inv(M)*m))

D=diag([0,lam,lam,lam,lam,lam,lam,lam,lam]);
A=inv(M)*D*M;
iD=diag([0,1/lam,1/lam,1/lam,1/lam,1/lam,1/lam,1/lam,1/lam]);
iA=inv(M)*iD*M;


clc
clear all
close all

nN = 20;
nB = 15;

MM(1) = 1.0d0;
MM(2) = 2.0d0;
MM(3) = 3.0d0;
phi = 1./MM;
a   = 0.0;

% User data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NU  = NU(nN);
B13 = B(MM(1),MM(3),nB);
B23 = B(MM(2),MM(3),nB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p1 = dlmread('p1-20-15.csv',',');
p2 = dlmread('p2-20-15.csv',',');
p3 = dlmread('p3-20-15.csv',',');
p  = p1+p2+p3;

r1 = dlmread('r1-20-15.csv',',');
r2 = dlmread('r2-20-15.csv',',');
r3 = dlmread('r3-20-15.csv',',');
r  = r1+r2+r3;

u1 = dlmread('u1-20-15.csv',',');
u2 = dlmread('u2-20-15.csv',',');
u3 = dlmread('u3-20-15.csv',',');

nx=size(p1,1);
ny=size(p1,2);
x=[1:nx];
y=[1:ny];

%figure(1)
%if (size(p1,1)==1)
%    plot(p1),grid on,hold on
%    plot(p2),grid on,hold on
%    plot(p3),grid on,hold on
%else
%    mesh(p1),grid on,hold on
%    mesh(p2),grid on,hold on
%    mesh(p3),grid on,hold on
%end

%xlabel('x axis')
%ylabel('y axis')
%zlabel('partial pressures')

y1 = p1./p;
y2 = p2./p;
y3 = p3./p;

x1 = r1./r;
x2 = r2./r;
x3 = r3./r;

v = y1.*u1+y2.*u2+y3.*u3;
u = x1.*u1+x2.*u2+x3.*u3;

iP = round(size(u1,1)/2);
jP = round(size(u1,2)/2);

gy1P = ( y1(iP,jP+1)-y1(iP,jP-1) )/2;
gy2P = ( y2(iP,jP+1)-y2(iP,jP-1) )/2;
gy3P = ( y3(iP,jP+1)-y3(iP,jP-1) )/2;

%figure(4)
%if (size(p1,1)==1)
%plot(y1(1,:).*u1(1,:),'b'),grid on,hold on
%plot(y2(1,:).*u2(1,:),'r'),grid on,hold on
%plot(y3(1,:).*u3(1,:),'k'),grid on,hold on

fh2 = figure(2);
load analytical

hl = plot(t,y(:,1),'k-');grid on,hold on
set(hl,'LineWidth',2,'MarkerSize',12)    
hl = plot(t,y(:,2),'k-.');grid on,hold on
set(hl,'LineWidth',2,'MarkerSize',12)    
hl = plot(t,1-y(:,1)-y(:,2),'k:');grid on,hold on
set(hl,'LineWidth',2,'MarkerSize',12)    

hl = plot(y1(1,:),'ko');hold on, grid on
set(hl,'LineWidth',2,'MarkerSize',12)    
hl = plot(y2(1,:),'ks');hold on, grid on
set(hl,'LineWidth',2,'MarkerSize',12)    
hl = plot(y3(1,:),'kv');hold on, grid on
set(hl,'LineWidth',2,'MarkerSize',12)    

axis([1 60 0 1])
hxl=xlabel('x / \delta x [-]');
hyl=ylabel('Molar Concentration y_{\sigma} [-]');
set(hyl,'FontSize',18,'FontWeight','bold')
set(hxl,'FontSize',18,'FontWeight','bold')
ax1 = gca;
set(ax1,'XColor','k','YColor','k',...
           'FontSize',18,...
           'LineWidth',1);
set(fh2,'Position',[360 41 728 673])
legend('Theory Species 1',...
    'Theory Species 2',...
    'Theory Species 3',...
    'Simulations Species 1',...
    'Simulations Species 2',...
    'Simulations Species 3')       
title('Stefan Tube N_3 = -6.1776 \cdot 10^{-5}')

%figure(3)
%plot(r1(1,:),'b'),hold on, grid on
%plot(r2(1,:),'r'),hold on, grid on
%plot(r3(1,:),'r'),hold on, grid on
%xlabel('x axis')
%ylabel('y axis')

%else
%    mesh(u1),grid on,hold on
%    mesh(u2),grid on,hold on
%    mesh(u3),grid on,hold on
%end

%figure(4)
%if (size(p1,1)==1)
%    subplot('121'),plot(u3),grid on,hold on
%    subplot('122'),plot(p3),grid on,hold on
%else
%    subplot('121'),mesh(u3),grid on,hold on
%    subplot('122'),mesh(p3),grid on,hold on
%end

%figure(5)
%plot([1:40]',u(:,30)),grid on,hold on

% Maxwell-Stefan model
%---------------------------------------
%Fomg1=phi(1)/(3*D(MM(1)))
%Fomg2=phi(2)/(3*D(MM(2)))
%Fomg3=phi(3)/(3*D(MM(3)))

%MSomg1=p(iP,jP)/r(iP,jP)*B(MM(1),MM(1))
%MSomg2=p(iP,jP)/r(iP,jP)*B(MM(2),MM(2))
%MSomg3=p(iP,jP)/r(iP,jP)*B(MM(3),MM(3))

% SOLVENT approximation
jM = size(p1,2);
gy1 = ( y1(iP,3:jM)-y1(iP,1:(jM-2)) )./2;
gy2 = ( y2(iP,3:jM)-y2(iP,1:(jM-2)) )./2;
gy3 = ( y3(iP,3:jM)-y3(iP,1:(jM-2)) )./2;

flux1 = y1(iP,2:(jM-1)).*(v(iP,2:(jM-1))-u1(iP,2:(jM-1)));
flux2 = y2(iP,2:(jM-1)).*(v(iP,2:(jM-1))-u2(iP,2:(jM-1)));

mmMSB13 = gy1(1:(jM-2))./( flux1 );
mmMSB23 = gy2(1:(jM-2))./( flux2 );

mmMSB13(find(mmMSB13<0))=0;
tent = mean(mmMSB13);
cj = 0;
mj = 0;
for j=1:size(mmMSB13,2);
    if abs(mmMSB13(j)-tent)/abs(tent)<1
        cj = cj+1;
        mj = mj+mmMSB13(j);
    end
end
%mMSB13 = mean(mmMSB13);
mMSB13 = mj/cj;

mmMSB23(find(mmMSB23<0))=0;
tent = mean(mmMSB23);
cj = 0;
mj = 0;
for j=1:size(mmMSB23,2);
    if abs(mmMSB23(j)-tent)/abs(tent)<1
        cj = cj+1;
        mj = mj+mmMSB23(j);
    end
end
%mMSB23 = mean(mmMSB23);
mMSB23 = mj/cj;

% Poiseuille flow
%---------------------------------------

Ny = size(u,1)+1;
uM = max(max(u));
uS = 2*u(1)-u(2);
mNU = a*Ny^2/(8*(uM-uS));

% Summary
%---------------------------------------

Sc1 = NU*B13;
Sc2 = NU*B23;

mSc1 = mNU*mMSB13;
mSc2 = mNU*mMSB23;

[nN,nB,NU,mNU,B13,mMSB13,B23,mMSB23,Sc1,mSc1,Sc2,mSc2]

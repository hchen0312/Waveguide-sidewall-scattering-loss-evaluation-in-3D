function [Green,detect] = G( rc,rcp,d1,d2,es,ea,lambda )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% The structure of RT_g is here: [R12,R23,T12,T23,T13]=RT_g(kza,kzs,ea,es,t,pol)

% calculate the stationary point at rp, rcp first.
xc=rc(1);
yc=rc(2);
zc=rc(3);
xcp=rcp(1);
ycp=rcp(2);
zcp=rcp(3);
dr=rc-rcp;
drl=sqrt(dr(1)^2+dr(2)^2+dr(3)^2);
k1=2*pi*sqrt(es)/lambda; % the wavenumber in sapphire
k3=k1;
k2=2*pi*sqrt(ea)/lambda; % wavenumber in aln
kx=k1*(xc-xcp)/drl; %stationary point of kx
ky=k1*(yc-ycp)/drl; %stationary point of ky
k1z=sqrt(k1*k1-kx*kx-ky*ky);
k2z=sqrt(k2*k2-kx*kx-ky*ky);
k3z=k1z;
% caculate the R and T
RT_tempTE=linspace(0,0,5);
RT_tempTM=linspace(0,0,5);
[RT_tempTE(1),RT_tempTE(2),RT_tempTE(3),RT_tempTE(4),RT_tempTE(5)]=RT_g(k2z,k1z,ea,es,d2-d1,true);
[RT_tempTM(1),RT_tempTM(2),RT_tempTM(3),RT_tempTM(4),RT_tempTM(5)]=RT_g(k2z,k1z,ea,es,d2-d1,false);
R23TE=RT_tempTE(2);
R21TE=RT_tempTE(2);
T21TE=RT_tempTE(4);
T23TE=RT_tempTE(4);
R23TM=RT_tempTM(2);
R21TM=RT_tempTM(2);
T21TM=RT_tempTM(4);
T23TM=RT_tempTM(4);
% calculate the A, B, and M
MTE=1-R23TE*R21TE*exp(2j*k2z*(d2-d1));
MTM=1-R23TM*R21TM*exp(2j*k2z*(d2-d1));

ATE=(exp(1j*(k1z-k2z)*zcp)+exp(1j*k2z*(zcp+2*d2)+1j*k1z*zcp)*R23TE)*T21TE;
ATM=(exp(1j*(k1z-k2z)*zcp)+exp(1j*k2z*(zcp+2*d2)+1j*k1z*zcp)*R23TM)*T21TM;

BTE=(exp(-1j*(k3z-k2z)*zcp)+exp(-1j*k2z*(zcp+2*d1)-1j*k3z*zcp)*R21TE)*T23TE;
BTM=(exp(-1j*(k3z-k2z)*zcp)+exp(-1j*k2z*(zcp+2*d1)-1j*k3z*zcp)*R21TM)*T23TM;

% construct CM and CNp, CNm
ks=[kx,ky,0];
k1p=[kx,ky,k1z];
k1m=[kx,ky,-k1z];
k2p=[kx,ky,k2z];
k2m=[kx,ky,k2z];
CM=cross(ks,[0,0,1])'*cross(ks,[0,0,1]);
CNp=cross(k1p,cross(ks,[0,0,1]))'*cross(k2m,cross(ks,[0,0,1]))/(k1*k2);
CNm=cross(k1m,cross(ks,[0,0,1]))'*cross(k2p,cross(ks,[0,0,1]))/(k1*k2);
%CM=[1 1 1;1 1 1;1 1 1];
%CNp=CM;
%CNm=CM;
% return the Green's function
if zc>0
    Green=exp(1j*dot(k1p,rc-rcp))*(k1z^2)*(CM*ATE/MTE+CNp*ATM/MTM)*exp(1j*(k1z-k2z)*d1)/(4*pi*(zc-zcp)*k1*norm(ks)*norm(ks)*k2z);
    %detect=exp(1j*k1z*(zc-zcp))*(ATE/MTE+ATM/MTM)*exp(1j*(k1z-k2z)*d1)/(4*pi*(zc-zcp)*k1*norm(ks)*norm(ks)*k2z);
    Temp=CM*[0 0 1]';
    detect=Temp(3);
elseif zc<-d2
    Green=-exp(1j*dot(k1m,rc-rcp))*(k3z^2)*(CM*BTE/MTE+CNm*BTM/MTM)*exp(-1j*(k3z-k2z)*d2)/(4*pi*(zc-zcp)*k3*norm(ks)*norm(ks)*k2z);
    %detect=exp(1j*k1z*(zcp-zc))*(BTE/MTE+BTM/MTM)*exp(-1j*(k3z-k2z)*d2)/(4*pi*(zc-zcp)*k3*norm(ks)*norm(ks)*k2z);
    Temp=CM*[0 0 1]';
    detect=Temp(3);
else
    Green=[0 0 0; 0 0 0; 0 0 0];
    detect=0;
end

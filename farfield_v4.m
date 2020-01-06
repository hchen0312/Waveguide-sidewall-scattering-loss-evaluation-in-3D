function [ E,S,Sr,power_p,power_sr,ratio,Er,P ] = farfield_v4(mm,Ex_mode,Ey_mode,Ez_mode,sl,dx,dy,h2,side,delta_epsilon,Green_label,dS,np,r,d1,d2,es,ea,lambda,omega,u,pol,sigma,Lc,beta)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% This function calculate the farfield poynting vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It needs:
% ii, figure number
% Ex,Ey,Ez_mode, the field distribution of mode
% sl, source location,
% dx,dy, resolution of mesh
% h2, waveguide height
% side, distance from boundry to waveguide
% Green_lable, ture for dislocation, false for sidewall
% dS, the 'dS' in E field integration
% np, points of far field
% r, distance from (0,0) to far field
% d1, d2, describe the slab thickness
% es, ea, epsilon of sapphire and aln
% lambda,omega, wavelength and angular frequency
% u, permibeality
% pol, simplified boundry, false indicate that the upper and lower boundry are neglected. 
% sigma, roughness in nm
% Lc, correlation length
% beta, propagation constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It will return
% E, far field E component
% S, far field poynting vector
% Sr, far field poynying radial component
% Power_p, total integration of poynting vector
% power_sr, total integration of radial component of poynting vector
% ratio, Power_p/Power_sr, should be equal to 1
% Er, radial component of E field, should be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It will plot
% figure ii: J distribution
% figure ii+1: far field dependence on theta
% figure ii+2: far field dependence on phi
% figure ii+3: surface plot of Sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Green_label
    index_x=round(side/dx+(d2*10^(6)/2/dx-sl*10^(6)/dx));
else
    index_x=round(side/dx);
end
index_yl=round(side/dy+1);
index_yu=round(side/dy+h2/dy);
cons=-1j*omega*delta_epsilon*8.854*10^(-12);
Jx=cons*Ex_mode(index_x,index_yl:index_yu);
%Jx=0;
Jy=cons*Ey_mode(index_x,index_yl:index_yu);
%Jy=0;
Jz=cons*Ez_mode(index_x,index_yl:index_yu);
%Jz=0;
figure(mm);
subplot(131);
plot(abs(Jx));
title('Jx');
subplot(132);
plot(abs(Jy));
title('Jy');
subplot(133);
plot(abs(Jz));
title('Jz');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=zeros(index_yu-index_yl+1,3); % This define the line current source;
ns=index_yu-index_yl+1; % points that describe the source
if Green_label % if the source is from dislocation
    if pol
        J(:,1)=Jx;
        J(:,3)=Jy;
        J(:,2)=Jz;
        rcpx=linspace(-d2/2-sl,-d2/2-sl,ns);
        rcpy=linspace(0,0,ns);
        rcpz=linspace(-h2*10^(-6),0,ns);
        %rcpz=linspace(-h2*10^(-6)/4,0,h2/dy);
        d2=h2*10^(-6);
    else 
        J(:,3)=Jx;
        J(:,2)=Jy;
        J(:,1)=Jz;
        rcpx=linspace(0,0,ns);
        rcpy=linspace(-h2*10^(-6),0,ns);
        rcpz=linspace(-d2/2-sl,-d2/2-sl,ns);
    end
else
    rcpz=linspace(0,0,ns);
    J(:,3)=Jx;
    J(:,2)=Jy;
    J(:,1)=Jz;
    rcpx=linspace(0,0,ns);
    rcpy=linspace(-h2*10^(-6),0,ns);
end
theta=linspace(0,pi,np);
phi=linspace(-pi,pi,np);
E=zeros(np,np,3);
rc=zeros(np,np,3);
detector=zeros(np,np);
S=zeros(np,np,3); % far field poynting vector
for ii=1:np
    for jj=1:np
        rc(ii,jj,1)=r*sin(theta(ii))*cos(phi(jj));
        rc(ii,jj,2)=r*sin(theta(ii))*sin(phi(jj));
        rc(ii,jj,3)=r*cos(theta(ii));
    end
end

for ii=1:np
    for jj=1:np
        for kk=1:ns
            rx=rc(ii,jj,1);
            ry=rc(ii,jj,2);
            rz=rc(ii,jj,3);
            if Green_label
                [Gtemp,detector(ii,jj)]=G([rx,ry,rz],[rcpx(kk),rcpy(kk),rcpz(kk)],d1,d2,es,ea,lambda);
            else
                [Gtemp,detector(ii,jj)]=G_sidewall([rx,ry,rz],[rcpx(kk),rcpy(kk),rcpz(kk)],d1,d2,es,ea,lambda);
            end
            Jtemp=J(kk,:);
            Etemp1=E(ii,jj,:);
            Etemp2=1j*omega*u*Gtemp*Jtemp';
            E(ii,jj,1)=Etemp1(1)+Etemp2(1);
            E(ii,jj,2)=Etemp1(2)+Etemp2(2);
            E(ii,jj,3)=Etemp1(3)+Etemp2(3);
            % function Green = G( rc,rcp,d1,d2,es,ea,lambda )
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% post treatment goes here
E=E/ns*h2*10^(-6)*dS;
Ex=E(:,:,1);
Ey=E(:,:,2);
Ez=E(:,:,3);
Er=zeros(np,np);
Sr=zeros(np,np); %at far field, E and H fields are fully tangential, and the Sr should euqal to total S
itg=zeros(np,np);
Ephi=zeros(np,np);
Etheta=zeros(np,np);
%Compute the far field poynting vector
k1=2*pi*sqrt(es)/lambda; % wavenumber in sapphire
for ii=1:np
    for jj=1:np
        S(ii,jj,:)=0.5*real(cross(E(ii,jj,:),cross(k1/(omega*u)*rc(ii,jj,:)/sqrt(sum(rc(ii,jj,:).*rc(ii,jj,:))),conj(E(ii,jj,:)))));
        Sx=S(ii,jj,1);
        Sy=S(ii,jj,2);
        Sz=S(ii,jj,3);
        Er(ii,jj)=Ez(ii,jj)*cos(theta(ii))+Ex(ii,jj)*cos(phi(jj))*sin(theta(ii))+Ey(ii,jj)*sin(phi(jj))*sin(theta(ii));
        Etheta(ii,jj)=-Ez(ii,jj)*sin(theta(ii))+Ex(ii,jj)*cos(theta(ii))*cos(phi(jj))+Ey(ii,jj)*cos(theta(ii))*sin(phi(jj));
        Ephi(ii,jj)=-sin(phi(jj))*Ex(ii,jj)+cos(phi(jj))*Ey(ii,jj);
        Sr(ii,jj)=Sz*cos(theta(ii))+Sx*cos(phi(jj))*sin(theta(ii))+Sy*sin(phi(jj))*sin(theta(ii));
        itg(ii,jj)=Sr(ii,jj)*R(sigma,Lc,beta-sqrt(es)*2*pi/lambda*sin(theta(ii))*cos(phi(jj)));
        % we also calculate the poynting vector on r direction, to verify the far field property
    end
end
Ss=sqrt(S(:,:,1).^2+S(:,:,2).^2+S(:,:,3).^2); % This is the total sum of poynting vector
Ss(find(isnan(Ss)))=0;
Sr(find(isnan(Sr)))=0;
itg(find(isnan(itg)))=0;
d_theta=theta(2)-theta(1);
d_phi=phi(2)-phi(1);
dS=r*r*sin(theta)*d_theta*d_phi;
itg=itg';
Ss=Ss';
Sr=Sr';
power_p=sum(sum(Ss).*dS);
power_sr=sum(sum(Sr).*dS);
P=sum(sum(itg).*dS);
ratio=power_p/power_sr; % should be 100%, if poynting vector is fully propagating
%%%%%% plotting
Sr=Sr';
Sphi=Sr(np/2+1,:);
Sr=Sr';
Stheta=Sr(np/2+1,:);
%figure(mm+1)
%polar(phi,Sphi);
%figure(mm+2)
%polar(theta,Stheta);
if Green_label
else
figure(mm+1)
Sr=Sr';
[az,el]=meshgrid(phi+pi,theta-pi/2);
%[az,el]=meshgrid(phi,theta);
[x_p,y_p,z_p]=sph2cart(az,el,Sr);
surf(x_p,y_p,z_p);
%surf(abs(Sr));
%view([-90 90])
end
end
% Please read "README.md" before runing this code

% The example computes scattering loss from TE00 mode within AlN waveguide
% WG geometry: H = 0.5 um, W = 1.2 um
clear all

nd2=1; % Number of width to iterate, setting to 1 to disable the loop
nh=1; % Number of height to iterate, setting to 1 to disable the loop
height=linspace(0.5,0.5,nh); % Waveguide height (in um)
d2record=linspace(1.2*10^(-6),1.2*10^(-6),nd2); % Waveguide width (in m)
dB_sidewall=zeros(nh,nd2);
%dB_dl_array=zeros(nh,nd2);
lambda=0.8*10^(-6); % wavelength defined in m
ea=2.1*2.1; % epsilon of aln
es=1.73*1.73; % epsilon of sapphire
sigma=1*10^(-9); % sidewall roughness, usually 3-5 nm
Lc=50*10^(-9); % correlation length, usually 50-100 nm
np=100; % number of points used to describe E field 
dd=1*10^(10); % dislocation density in cm-2

for ii=1:nh
for jj=1:nd2
as_ratio=d2record(jj)/(height(ii)*10^(-6));
d1=0;
d2=d2record(jj); % d2-d1 equal to the waveguide width, in um.
r=10000*lambda;
omega=2*pi*3*10^(8)/lambda;
u=4*pi*10^(-7);
%sl=0*10^(-6); %It defines the place of dislocation away from center of waveguide d/2

%sl=linspace(0,d2/2,11);
sl=d2/4;
%power_dl_rc=linspace(0,0,11);
%power_sd_rc=linspace(0,0,11);
%Sr_dl_rcd=zeros(100,100,11);
delta_epsilon_dl=1; %It defines the pertubation of epsilon
delta_epsilon_sd=ea-es;
%delta_epsilon_sd=2;
dS_dl=1*10^(-18);
dS_sd=dS_dl;
%dS_sd=1*10^(-18);
%dS_dl=1*10^(-18);
% To define current density, the mode is calculated in this section
%%%%%%%%%%%%%%%%%%%% mode calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmodes=12; % TE and TM
TE=1;
if as_ratio>=1
    if TE==1
        order=1;
    else
        order=2;
    end
else
    if TE==1
        order=2;
    else
        order=1;
    end
end
%order=1 % NOTE: we can get certain mode by either setting "order" or symmetry from "000S" to "000A"
pole=false; % This term is related to the simplification of modal. In most of the situation, it should be false.
% 1st order, TE
% 2nd order, TM
% Refractive indices:
n1=1.73;
n2= 2.1;
n3 = n1;          % Upper cladding (alloy of SiO2 and SiN)
% Layer heights:
h2=height(ii);
%h2 =height(ii);           % Core thickness
h3 = 1;           % Upper cladding
% Horizontal dimensions:
rh = h2;           % Ridge height
%rw=1.5;
rw = d2*10^(6)/2;           % Ridge half-width
%THIS IS RELATED TO d2!
side = 1;         % Space on side
% Grid size:
dx = 0.01;        % grid size (horizontal)
dy = 0.005;        % grid size (vertical)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,xc,yc,nx,ny,eps,edges] = ...
waveguidemeshfull([n1,n2,n3],[side,h2,h3],rh,rw, ...
                  side,dx,dy);
[Hx_mode,Hy_mode,ntemp] = wgmodes(lambda*10^(6),n2,nmodes,dx,dy,eps,'0000');
neff=ntemp(order);
beta=neff*2*pi/lambda;
[Hz_mode,Ex_mode,Ey_mode,Ez_mode] = postprocess (lambda*10^(6), neff, Hx_mode(:,:,order), Hy_mode(:,:,order), dx, dy, eps, '0000');
fprintf(1,'neff = %.6f\n',neff);
%[Ex_mode,Ey_mode,Ez_mode,Hx_mode,Hy_mode,Hz_mode] = normalize(dx,dy,Ex_mode,Ey_mode,Ez_mode,Hx_mode,Hy_mode,Hz_mode);
figure(1);
subplot(131);
contourmode(x,y,Ex_mode(:,:));
title('Ex');
for v = edges, line(v{:}); end
subplot(132);
contourmode(x,y,Ey_mode(:,:));
title('Ey');
for v = edges, line(v{:}); end
subplot(133);
contourmode(x,y,Ez_mode(:,:));
title('Ez');
for v = edges, line(v{:}); end

%%%comput the power of mode
[nx_mode,ny_mode] = size(Ex_mode);
S_mode=zeros(nx_mode,ny_mode);
Z=377/2;
for iiii=1:nx_mode
    for jjjj=1:ny_mode
        if TE==1
            E=Ex_mode(iiii,jjjj);
        else
            E=Ey_mode(iiii,jjjj);
        end
        S_temp=0.5*real(E*conj(E/Z));
        S_mode(iiii,jjjj)=S_temp;
    end
end
P=sum(sum(S_mode))*dx*dy*10^(-12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%[E_dl,S_dl,Sr_dl_rcd,power_p_dl,power_sr_dl,ratio_dl,Er_dl,P_dB_dl ]=farfield_v4(2,Ex_mode,Ey_mode,Ez_mode,sl,dx,dy,h2,side,delta_epsilon_dl,true,dS_dl,np,r,d1,d2,es,ea,lambda,omega,u,pole,sigma,Lc,beta);
%disp(['power_p_dl=',num2str(power_p_dl)]);
%disp(['power_sr_dl=',num2str(power_sr_dl)]);
[E_sd,S_sd,Sr_sd,power_p_sd,power_sr_sd,ratio_sd,Er_sd,P_dB ]=farfield_v4(2,Ex_mode,Ey_mode,Ez_mode,sl,dx,dy,h2,side,delta_epsilon_sd,false,dS_sd,np,r,d1,d2,es,ea,lambda,omega,u,pole,sigma,Lc,beta);
disp(['power_p_sd=',num2str(power_p_sd)]);
disp(['power_sr_sd=',num2str(power_sr_sd)]);
%Efficiency_ratio=power_p_dl/power_p_sd;
P_dB=P_dB*10^(27);
%P_dB_dl=P_dB_dl*10^(27);
%power_dl_rc=power_p_dl;
power_sd_rc=power_p_sd;
%power_dl_rc=power_dl_rc';
power_sd_rc=power_sd_rc';
%ratio=power_dl_rc./power_sd_rc;
dB_sidewall(ii,jj)=2*10*log10((P-P_dB)/P)*10^(7);
%dB_dl_array(ii,jj)=10*log10((P-P_dB_dl)/P)*10^(7);
end
end
dB_sidewall_real=dB_sidewall*(ea-es);
disp(['loss=',num2str(dB_sidewall_real),'dB/cm']);

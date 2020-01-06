function [R12,R23,T12,T23,T13]=RT_g(kza,kzs,ea,es,t,pol)
% kza, wavevector along z direction in aln layer
% kzs, wavevector along z direction in sapphire layer
% ea, eps of aln layer
% es, eps of sapphire layer
% this function returns generalized reflection and transmission at 12 and 23 interface, by
% calling R_TM and R_TE
% t, thickness (width of waveguide)
% TE: pol=ture, TM: pol=false
if (pol)
    R12=(R_TE(kzs,kza)+R_TE(kza,kzs).*exp(2*1j*kza*t))./(1+R_TE(kzs,kza).*R_TE(kza,kzs).*exp(2*1j*kza*t));
    R23=R_TE(kza,kzs);
    T12=1+R12;
    T23=1+R23;
    T13=(T12./(1+R_TE(kzs,kza).*R23.*exp(2*1j*kza*t))).*(T23);
else
    R12=(R_TM(kzs,kza,es,ea)+R_TM(kza,kzs,ea,es).*exp(2*1j*kza*t))./(1+R_TM(kzs,kza,es,ea).*R_TM(kza,kzs,ea,es).*exp(2*1j*kza*t));
    R23=R_TM(kza,kzs,ea,es);
    T12=1+R12;
    T23=1+R23;
    T13=(T12./(1+R_TM(kzs,kza,ea,es).*R23.*exp(2*1j*kza*t))).*(T23);
end
end


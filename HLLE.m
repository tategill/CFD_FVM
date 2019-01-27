function [F_hat,Se,error] = HLLE(uL, uR, n, gamma)
%HLLE computes the numerical HLLE flux for states uL and uR where the unit
%normal n (row vector) points from left to right, and gamma is the ratio of
%specific heats.
%   Addtionally HLLE computes the maximum wave speed S.

error = 0; %by deafult

[FLx,FLy,VL,cL,errorL] = eulerflux(uL, gamma);
[FRx,FRy,VR,cR,errorR] = eulerflux(uR, gamma);

if errorL == 1 || errorR == 1 %If recieve error from eulerflux, pass error on 
    error = 1;
end
    
vL = dot(VL,n); vR = dot(VR,n);
FL = FLx*n(1) + FLy*n(2);
FR = FRx*n(1) + FRy*n(2);

sLmin = min([0 vL-cL]); sRmin = min([0 vR-cR]);
sLmax = max([0 vL+cL]); sRmax = max([0 vR+cR]);
smin = min([sLmin sRmin]); smax = max([sLmax sRmax]); 

F_hat = 0.5*(FR+FL) - 0.5*((smax+smin)/(smax-smin))*(FR-FL) + ((smax*smin)/(smax-smin))*(uR-uL);

Se = max([abs(vL)+cL abs(vR)+cR]);

end


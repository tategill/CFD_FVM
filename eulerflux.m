function [Fx,Fy,V,c,error] = eulerflux(u,gamma)
%eulerflux calculates the 2D flux components via the euler equations from 
%the state u and the ratio of specific heats gamma.
%   additional outputs are the velocity vector of the state V, and the
%   sound speed c.
%u = [rho, rho*u, rho*v, rho*E];
%Fx = [rho*u, rho*u^2+P, rho*u*v, rho*u*H];
%Fy = [rho*v, rho*v*u, rho*v^2+P, rho*v*H];
%V = [Vx Vy];
%c = srqt(gamma*P/rho);

error = 0; %by deafult

V = [u(2)/u(1) u(3)/u(1)];
magV2 = (norm(V))^2;
P = (gamma-1)*(u(4) - 0.5*u(1)*magV2);

if P < 0 %throw error if pressure is negative 
    error = 1;
end

H = u(4)/u(1) + P/u(1);
c = sqrt(gamma*P/u(1));

Fx = [u(2) u(2)^2/u(1) + P (u(2)*u(3))/u(1) u(2)*H];
Fy = [u(3) (u(2)*u(3))/u(1) u(3)^2/u(1) + P u(3)*H];

end


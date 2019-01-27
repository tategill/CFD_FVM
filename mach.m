function M = mach(u, gamma)
%mach calclates the mach number M for the euler state u and the ratio of 
%specific heats gamma

velocity = 1./u(:,1).*sqrt(u(:,2).^2+u(:,3).^2); %magnitude of velocity
pressure = (gamma-1).*(u(:,4)-0.5.*u(:,1).*velocity.^2); 
soundspeed = (gamma*pressure./u(:,1)).^(0.5);
M = velocity./soundspeed; %mach number
end


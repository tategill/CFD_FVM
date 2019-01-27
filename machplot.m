function machplot(u,E,V,gamma,heading,zoom)
%machplot creates a 2D color plot of the local mach numbers for
%the state matrix u.
%   gamma is thre ratio of specific heats, used in the mach calculation
%   heading is a string that becomes the title of the figure
%   zoom is the zoom factor (1x 5x 10x ect) about the origin

%CALCULATE MACH NUMBER
velocity = 1./u(:,1).*sqrt(u(:,2).^2+u(:,3).^2); %magnitude of velocity
pressure = (gamma-1).*(u(:,4)-0.5.*u(:,1).*velocity.^2); 
soundspeed = (gamma*pressure./u(:,1)).^(0.5);
mach = velocity./soundspeed; %local mach number at every point

meshplot(mach,E,V,'Local Mach Number',heading,zoom,0);

end




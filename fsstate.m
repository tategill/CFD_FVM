function ufs = fsstate(M,alpha,gamma)
%fsstate calculates the free stream state given a mach number M, an angle
%of attack alpha in degrees, and a ratio of specific heats gamma

ufs = zeros(1,4);

ufs(1) = 1;
ufs(2) = M*cosd(alpha);
ufs(3) = M*sind(alpha);
ufs(4) = 1/((gamma-1)*gamma) + (M^2)/(2);

end
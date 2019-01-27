function [bool] = isareapos(V,p1,p2,p3)
%isareapos calculates the area of a triangle from points 1 2 and 3, and
%returns 1 is the area is positive and zero if the area is negative.
%   points (p1,p2,p3) are [x y] coordinates of the triangle's verticies.

x = [V(p1,1); V(p2,1); V(p3,1)];
y = [V(p1,2); V(p2,2); V(p3,2)];          
area = (x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)))/2;

if area >= 0     
    bool = 1;
elseif area < 0
    bool = 0;
end  

end
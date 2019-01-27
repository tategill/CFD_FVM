function meshplot(u,E,V,colorlabel,heading,zoom,nodat)
%Generates contour plot of data in u onto mesh defined by E and V on
%current plot
%   u must be a vector with length = to length of E
%   colorlabel is a string that is assigned to the colorbar scale.
%   heading is a string that acts as the title of the plot (could be any
%   other desired string)
%   zoom is the zoom factor (1x 5x 10x ect) about the origin
%   nodat is a boolean to turn either the lines of the mesh on or off and to plot white tirangles. 

hold on

set(gca,'FontSize',12)

%Set initial axis to the range of V and be unstreched
ax0 = [min(V(:,1)) max(V(:,1)) min(V(:,2)) max(V(:,2))]; 
axis(ax0,'equal');

%Set zoomed axis to be reduced by zoom factor
ax = axis./(zoom);
axis(ax);

cbar = colorbar;
cbar.Label.String = colorlabel;
xlabel('x Position'); ylabel('y Position'); 
title(heading,'Interpreter','none');

%Get approx plotting limits for zoom axis
limx_neg = ax(1)*1.75;
limx_pos = ax(2)*1.75;
limy_neg = ax(3)*1.75;
limy_pos = ax(4)*1.75;

%DRAW TRIANGLES
for i = (1:length(E(:,1)))
    
    tri = zeros(3,2); %matrix that defines x and y cordinates for each point of triangle
    
    p1 = E(i,1); %points are the values of a particular row of E
    p2 = E(i,2);
    p3 = E(i,3); 
     
    tri(1,1) = V(p1,1); %assign x coordinates
    tri(2,1) = V(p2,1);
    tri(3,1) = V(p3,1);
 
    tri(1,2) = V(p1,2); %assign y coordinates
    tri(2,2) = V(p2,2);
    tri(3,2) = V(p3,2);
    
    %condition checks to make sure triagle is near the desired plot window
    if min(tri(:,1)) >= limx_neg && max(tri(:,1)) <= limx_pos && min(tri(:,2)) >= limy_neg && max(tri(:,2)) <= limy_pos
        if nodat
            fill(tri(:,1),tri(:,2),[1 1 1]);
        else
            fill(tri(:,1),tri(:,2),u(i),'LineStyle','none'); 
        end
    end
end

hold off

end
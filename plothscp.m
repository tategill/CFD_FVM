function plothscp(V,B,Cp,color)
%plothscp plots (on the current figure) the pressure coefficient as a
%function of the heat shield angle, as defined in the project description.
%   V is the vertex matrix for the mesh
%   B is the boundary info cell array
%   Cp is the converged Cp values on the capsule
%   color is the color number for the plot

%extract capsule boundary data
capB = B{5,3};

%get points corrensponding to capsule edges
x = [V(capB(:,1),1) V(capB(:,2),1)];
y = [V(capB(:,1),2) V(capB(:,2),2)];

%pair Cp with geometery data before sorting
pts_Cp = [Cp x(:,1) y(:,1) x(:,2) y(:,2)];

%just get edges on heat shield (i.e both x coords of edge points <=0) and
%sort bottom to top.
pts_hs = sortrows(pts_Cp(pts_Cp(:,2)<=0 & pts_Cp(:,4)<=0,:),3); 

%pts_hs = pts_Cp %Code for computing full Cp plots

circ_ctr = [.8 0]; %x and y location of heat sheild circle center

%calculate the rise and run from the circle center for starting(1) and
%ending points(2)
rise1 = (pts_hs(:,3)-circ_ctr(2));
run1 = -(pts_hs(:,2)-circ_ctr(1));
hyp1 = sqrt(rise1.^2 + run1.^2);

rise2 = (pts_hs(:,5)-circ_ctr(2));
run2 = -(pts_hs(:,4)-circ_ctr(1));
hyp2 = sqrt(rise2.^2 + run2.^2); 

%calculate angle and add to info matrix
for i = 1:size(pts_hs,1)
    if rise1(i) >= 0 
    pts_hs(i,6) = acosd(run1(i)/hyp1(i));  

    elseif rise1(i) <0 
    pts_hs(i,6) = -acosd(run1(i)/hyp1(i));
    
    end
    
    if rise2(i) >=0
    pts_hs(i,7) = acosd(run2(i)/hyp2(i));  
    
    elseif rise2(i) <0
    pts_hs(i,7) = -acosd(run2(i)/hyp2(i)); 
    
    end
end


pts_hs = sortrows(pts_hs,6);


%set up plotting data
theta = pts_hs(1,6):0.1:pts_hs(size(pts_hs,1),7);
Cp_theta = zeros(size(theta));

%find the values of theta that are between the starting and ending angles,
%and set the value of Cp(theta) as the Cp on that edge.
for i = 1:size(pts_hs,1)
    for t = 1:size(theta,2)
        if theta(t) >= pts_hs(i,6) && theta(t) < pts_hs(i,7)
            Cp_theta(t) = pts_hs(i,1);
        end
    end
end

load('colour.mat'); %custom color pallet 

hold on
plot(theta,Cp_theta,'Color',colour(color,:),'LineWidth',2)
line([0 0], ylim,'color',[0.5 0.5 0.5],'linewidth',2) %y-axis
grid on

%% Commented out code for computing full Cp plots on capsule
%{
%pts_hs(size(pts_hs,1),:) = []; %delete back one to eliminate error 

%set up plotting data for top/bottom
theta_pos = 0:0.001:180;
Cp_theta_pos = zeros(size(theta_pos));

%find the values of theta that are between the starting and ending angles,
%and set the value of Cp(theta) as the Cp on that edge.
for i = 1:size(pts_hs,1)
    for t = 1:size(theta_pos,2)
        if theta_pos(t) >= pts_hs(i,6) && theta_pos(t) < pts_hs(i,7)
            Cp_theta_pos(t) = pts_hs(i,1);
        end
    end
end

theta_neg = 0:-0.001:-180;
Cp_theta_neg = zeros(size(theta_neg));

%find the values of theta that are between the starting and ending angles,
%and set the value of Cp(theta) as the Cp on that edge.
for i = 1:size(pts_hs,1)
    for t = 1:size(theta_neg,2)
        if theta_neg(t) >= pts_hs(i,6) && theta_neg(t) < pts_hs(i,7)
            Cp_theta_neg(t) = pts_hs(i,1);
        end
    end
end

%set up plotting data for front/back
theta_front = -90:0.001:90;
Cp_theta_front = zeros(size(theta_front));

%find the values of theta that are between the starting and ending angles,
%and set the value of Cp(theta) as the Cp on that edge.
for i = 1:size(pts_hs,1)
    for t = 1:size(theta_front,2)
        if theta_front(t) >= pts_hs(i,6) && theta_front(t) < pts_hs(i,7)
            Cp_theta_front(t) = pts_hs(i,1);
        end
    end
end

theta_back1 = 90:+0.001:180;
theta_back2 = -180:+0.001:-90;
Cp_theta_back1 = zeros(size(theta_back1));
Cp_theta_back2 = zeros(size(theta_back2));

%find the values of theta that are between the starting and ending angles,
%and set the value of Cp(theta) as the Cp on that edge.
for i = 1:size(pts_hs,1)
    for t = 1:size(theta_back1,2)
        if theta_back1(t) >= pts_hs(i,6) && theta_back1(t) < pts_hs(i,7)
            Cp_theta_back1(t) = pts_hs(i,1);
        end
    end
    for t = 1:size(theta_back2,2)
        if theta_back2(t) >= pts_hs(i,6) && theta_back2(t) < pts_hs(i,7)
            Cp_theta_back2(t) = pts_hs(i,1);
        end
    end
end


%plot data
figure(1)
clf(1)

hold on
grid on

plot(theta_pos,Cp_theta_pos,'b','LineWidth',2)
plot(-theta_neg,Cp_theta_neg,'r','LineWidth',2)

line([0 0], ylim,'color','k','linewidth',2) %y-axis
xlabel('Angle from HS center [deg]'); ylabel('Pressure Coefficient Cp');

legend('Capsule Top','Capsule Bottom');

figure(2)
clf(2)

hold on
grid on

plot(theta_front,Cp_theta_front,'b','LineWidth',2)
plot(-theta_back1+180,Cp_theta_back1,'r',-theta_back2-180,Cp_theta_back2,'r','LineWidth',2)

line([0 0], ylim,'color','k','linewidth',2) %y-axis
xlabel('Angle from HS center [deg]'); ylabel('Pressure Coefficient Cp');

legend('Capsule Front','Capsule Back');
%}
end


%Visualize boundaries and normal directions

% Test Data 1
%{
E = [1,2,4;2,5,4;2,3,5;3,6,5];
V = [0,0;1,0;2,0;0,1;1,1;2,1];
B = cell(4,3);
B{1,1} = 'Left';    B{1,2} = 1; B{1,3} = [4 1];
B{2,1} = 'Top';     B{2,2} = 2; B{2,3} = [6 5; 5 4];
B{3,1} = 'Right';   B{3,2} = 1; B{3,3} = [3 6];
B{4,1} = 'Bottom';  B{4,2} = 2; B{4,3} = [1 2; 2 3];
%}

load('colour.mat');

gamma = 1.3;
[V, E, B] = readgri('Mesh0.gri');

%get boundary egde data
bedges = cell(size(B,1),1);
for j = 1:size(B,1)
    bedges{j} = bedgedat(E,V,B{j,3});
end

%get interior egde data
C = connect(E);
inedges = inedgedat(E,V,C);

%% Normal Directions
figure(1)
clf
hold on
axis equal
for j = 1:4
    for i = 1:size(B{j,3},1)
        line = [V(B{j,3}(i,1),:);V(B{j,3}(i,2),:)];
        plot(line(:,1),line(:,2),'Color',colour(3,:),'Linewidth',2)  
        mp = [mean(line(:,1));mean(line(:,2))];
        quiver(mp(1),mp(2),bedges{j}(i,3),bedges{j}(i,4),'Color',[0 0 0],'LineWidth',1);
    end
end
xlabel('x Position'); ylabel('y Position'); 
title('Normals for Exterior Boundaries on Initial Mesh')
set(gca,'FontSize',12)

figure(2)
clf
hold on
axis equal
for i = 1:size(B{5,3},1)
    line = [V(B{5,3}(i,1),:);V(B{5,3}(i,2),:)];
    plot(line(:,1),line(:,2),'Color',colour(3,:),'Linewidth',2)  
    mp = [mean(line(:,1));mean(line(:,2))];
    quiver(mp(1),mp(2),0.25*bedges{5}(i,3),0.25*bedges{5}(i,4),'Color',[0 0 0],'LineWidth',1);
end
xlabel('x Position'); ylabel('y Position'); 
title('Normals for Capsule Boundaries on Initial Mesh (1/4 length)')
set(gca,'FontSize',12)

%% HLLE TEST CASE

%Test Data 2 (HLLE FLUX)
E = [E(1,:);E(81,:)];

cntr = [0 0; 0 0]; %centroids of L and right elements

figure(3)
clf
hold on
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

    fill(tri(:,1),tri(:,2),colour(i+1,:),'LineStyle','-');
    
    cntr(i,1) = mean(tri(:,1));
    cntr(i,2) = mean(tri(:,2));
end



%First Edge [288, 15, 0.994631976904759, -0.103475748456005, 3.05529135838450, 1,  81, 0]
%           [nA,  nB, nx,                ny,                 dl,               EL, ER, flag]

line = [V(inedges(1,1),:);V(inedges(1,2),:)];
mp = [mean(line(:,1));mean(line(:,2))];
quiver(mp(1),mp(2),inedges(1,3),inedges(1,4),'k','LineWidth',2);


text(cntr(1,1),cntr(1,2)-.25,'L','FontSize',24)
text(cntr(2,1),cntr(2,2)-.25,'R','FontSize',24)
xlabel('x Position'); ylabel('y Position'); 
title('Flux Setup for 1st Interior Edge');
set(gca,'FontSize',12)

close all

%% Test Case 1a: Subsonic Regular Normals
uL = fsstate(0.5,0,gamma);
uR = fsstate(0.55,3,gamma);

n = [inedges(1,3) inedges(1,4)];
[F_hat,S] = HLLE(uL, uR, n, gamma);

[Fx,Fy,V,c] = eulerflux(uL,gamma);
FeulerL = Fx*n(1) + Fy*n(2);
SeulerL = c + abs(dot(V,n));

[Fx,Fy,V,c] = eulerflux(uR,gamma);
FeulerR = Fx*n(1) + Fy*n(2);
SeulerR = c + abs(dot(V,n));

%% Test Case 1b: Subsonic Flipped Normals
uL = fsstate(0.5,0,gamma);
uR = fsstate(0.55,3,gamma);

n = -1.*[inedges(1,3) inedges(1,4)];
[F_hat,S] = HLLE(uL, uR, n, gamma);

[Fx,Fy,V,c] = eulerflux(uL,gamma);
FeulerL = Fx*n(1) + Fy*n(2);
SeulerL = c + abs(dot(V,n));

[Fx,Fy,V,c] = eulerflux(uR,gamma);
FeulerR = Fx*n(1) + Fy*n(2);
SeulerR = c + abs(dot(V,n));


%% Test Case 2a: Supersonic Regular Normals

uL = fsstate(2.5,0,gamma);
uR = fsstate(2.55,3,gamma);

n = [inedges(1,3) inedges(1,4)];
[F_hat,S] = HLLE(uL, uR, n, gamma);

[Fx,Fy,V,c] = eulerflux(uL,gamma);
FeulerL = Fx*n(1) + Fy*n(2);
SeulerL = c + abs(dot(V,n));

[Fx,Fy,V,c] = eulerflux(uR,gamma);
FeulerR = Fx*n(1) + Fy*n(2);
SeulerR = c + abs(dot(V,n));

%% Test Case 2b: Supersonic Flipped Normals

uL = fsstate(2.5,0,gamma);
uR = fsstate(2.55,3,gamma);

n = -1*[inedges(1,3) inedges(1,4)];
[F_hat,S] = HLLE(uL, uR, n, gamma);

[Fx,Fy,V,c] = eulerflux(uL,gamma);
FeulerL = Fx*n(1) + Fy*n(2);
SeulerL = c + abs(dot(V,n));

[Fx,Fy,V,c] = eulerflux(uR,gamma);
FeulerR = Fx*n(1) + Fy*n(2);
SeulerR = c + abs(dot(V,n));

%{
uL = fsstate(2.5,0,gamma);

R = 0;
n = [inedges(1,3) inedges(1,4)];
[F_hat,S] = HLLE(u, u, n, gamma);
R = R + F_hat*inedges(1,5);

[Fx,Fy,V,c] = eulerflux(u,gamma);
Feuler = Fx*n(1) + Fy*n(2);
Seuler = c + abs(dot(V,n));
%}

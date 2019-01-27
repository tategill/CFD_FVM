%% SET UP

load('colour.mat');
fname = input('Enter mesh file name: ','s');

[V, E, B] = readgri(fname);

%gamma = 1.3; %ratio of specific heats
%M0 = 1.3; %inital state mach number
%Mfs = 8; %free stream mach number

alphafs = input('Enter free stream angle of attack (deg): '); %free stream angle of attack in degrees

answer = input('Would you like to use a state restart file? (y/n): ','s');

if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes')) 
restartU = input('Enter restart state file name: ','s');
u0 = fscanf(fopen(restartU, 'r'),'%f',[4,size(E,1)])';
else 
u0 = ones(size(E,1),4).*(fsstate(1.3,alphafs,1.3)); %inital startup state
end

fprintf('\n<<<<< Run Start >>>>>\n');

%% FVM CODE
[u, N, Rnorm, Cl, Cd, Cm, Cp, e_ind]= FVM(fname,2000,u0,8,alphafs,1.3);
if N == 2000
fprintf('/n...WARNING: SOLUTION MAY NOT BE CONVERGED!...\n\n');
end

fprintf('<<<<< FVM Complete >>>>>\n\n');

%% Plots
answer = input('Do you want to generate solution plots? (y/n): ','s');
if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes')) 
saveprefix = input('Enter save file prefix for this run: ','s');
   
%% RESIDUAL PLOT
figure(1)
clf(1)
semilogy(Rnorm,'Color',colour(3,:),'LineWidth',2)
line(xlim, [10^-5 10^-5],'color','k','linewidth',1,'LineStyle','--') %convergence tolerance line 
text(0,8.5^-5,'Convergence Tolerance = 10^{-5}');
xlabel('Time Step Iterations'); ylabel('L_2 Norm of State Residuals');
title('L_2 Norm Convergence');
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(gcf,strcat(saveprefix,'ResidualConv.png'))
input('\nNext Plot...?');
close(1)
fprintf('--------> Residual Plot: Done!\n');

%% AERO COEFFICIENT PLOTS
f = figure(2);
clf(2)
f.OuterPosition = [20 278 1250 420];
subplot('Position',[.05 .125 0.275 .8]);
plot(Cl,'Color',colour(3,:),'LineWidth',2)
xlabel('Time Step Iterations'); ylabel('Lift Coefficient');
title('Lift Coefficient Convergence');

subplot('Position',[.375 .125 0.275 .8]);
plot(Cd,'Color',colour(3,:),'LineWidth',2)
xlabel('Time Step Iterations'); ylabel('Drag Coefficient');
title('Drag Coefficient Convergence');

subplot('Position',[.7 .125 0.275 .8]);
plot(Cm,'Color',colour(3,:),'LineWidth',2)
xlabel('Time Step Iterations'); ylabel('Moment Coefficient');
title('Moment Coefficient Convergence');

set(findall(gcf,'-property','FontSize'),'FontSize',12);
saveas(gcf,strcat(saveprefix,'AeroCoeffConv.png'))

input('Next Plot...?');
close(2)

fprintf('--------> Aero Coefficient Plots: Done!\n');

%% MACH PLOTS
f = figure(3);
clf(3)
%f.OuterPosition = [100 278 560*2 420]; %[x y width height]
%default position = [360 278 560 420] 

%{
subplot('Position',[.05 .1 0.25 .85]); %[left bottom width height]
machplot(u,E,V,1.3,'Converged Mach Field (1x Zoom)',1);
hold on 
plotcapsule(B,V);
%}

%subplot('Position',[.375 .1 0.6 .85]);
machplot(u,E,V,1.3,'Converged Mach Field (10x Zoom)',10);
hold on 
plotcapsule(B,V);

set(findall(gcf,'-property','FontSize'),'FontSize',12);
saveas(gcf,strcat(saveprefix,'ConvMach.png'))

input('Next Plot...?');
close(3)

fprintf('--------> Mach Plots: Done!\n');

%% Cp PLOT
figure(4)
clf(4)
plothscp(V,B,Cp,3)
xlabel('Heat Shield Angle, \theta [Deg]'); ylabel('Coefficient of Pressure');
title('Converged Cp on the Heat Shield');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
saveas(gcf,strcat(saveprefix,'CpPlot.png'))

input('Next Plot...?');
close(4)

fprintf('--------> Cp Plot: Done!\n\n');

end

%% Plot Mesh
answer = input('Would you like plot the mesh? (y/n): ','s');

if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes'))
figure(4)
clf(4)
hold on
nulldat = zeros(size(E,1));
meshplot(nulldat,E,V,'colorlabel','Full Mesh',1,1)
colorbar('off');
saveas(gcf,strcat(saveprefix,'MeshFull.png'))

figure(5)
clf(5)
hold on
meshplot(nulldat,E,V,'colorlabel','Zoomed Mesh',10,1)
colorbar('off');
saveas(gcf,strcat(saveprefix,'MeshZoom.png'))

input('Done...?');
close(4)
close(5)

end

answer = input('Would you like add converged aero coefficients to Datalog? (y/n): ','s');

if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes')) 
%% Add Converged Values to Data Log    
MeshN = input('Enter Mesh Number: ');
Datalog{MeshN+1,1} = MeshN;
Datalog{MeshN+1,2} = size(E,1);
Datalog{MeshN+1,3} = Cl(N-1);
Datalog{MeshN+1,4} = Cd(N-1);
Datalog{MeshN+1,5} = Cm(N-1);
Datalog{MeshN+1,6} = Cp;
end

answer = input('Would you like refine mesh? (y/n): ','s');

if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes')) 
%% Mesh Adapt
[V2,E2,B2,u2] = meshadapt(V,E,B,u,e_ind);

%saveas(gcf,strcat(saveprefix,'AdaptPlot.png'))
%close(10)

end

fprintf('\n<<<<< Run Complete >>>>>\n\n');

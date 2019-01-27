function [u, N, Rnorm, Cl, Cd, Cm, Cp, e_ind]= FVM(fname,Niter,u0,Mfs,alphafs,gamma)
%FVM implements a finite volume method for the given mesh data in fname,
%and state u0
%   Niter is the maximun number of iterations/time steps
%   Mfs is the freestream mach number
%   alphafs is the freestream angle of attack
%   gamma is the ratio of specific heats

%% MESH CONDITIONING

[V, E, B] = readgri(fname);

%get boundary egde data
bedges = cell(size(B,1),1);
for j = 1:size(B,1)
    bedges{j} = bedgedat(E,V,B{j,3});
end

%get interior egde data
C = connect(E);
inedges = inedgedat(E,V,C);

%% INITIALIZATION

CFLmax = .9; %Maximum CFL allowed for time stepping
ufs = fsstate(Mfs,alphafs,gamma); %free stream state

u = u0;
R = zeros(size(E,1),4); %residual vector
S = zeros(size(E,1),1); %wave speed vector
e_ind = zeros(size(inedges,1) + size(bedges{5},1),2); %error indicator for interior and capsule edges

Rnorm = zeros(Niter,1); %Log of residual norm 
Cl = zeros(Niter,1); %Log of lift coefficient
Cd = zeros(Niter,1); %Log of drag coefficient 
Cm = zeros(Niter,1); %Log of moment coefficient 

%Initialize Force and Moment on the Capsule
Fx = 0;
Fy = 0;
M = 0;

%% MAIN ITERATION LOOP
count_max = 100; %counter parameters for periodic state output
count = 0;
CFLTuner = 0; %initialize CFL Tuner in the off condition

for N = 1:Niter
    %% SINGLE TIME STEP 
    R = R.*0; %initialize residuals 
    S = S.*0; %initialize wave speeds
    e_ind = e_ind.*0; %initialize error indicator 
    error = 0; %initialize error toggle for unreal states
    
    count = count+1; %add value to counter
    
    %% interior edges
    for e = 1:size(inedges,1)        
        if error == 1 %if there was an error
            break
        end
            
        %extract data: [nA nB nx ny dl EL ER]
        EL = inedges(e,6);
        ER = inedges(e,7);
        n = [inedges(e,3) inedges(e,4)];
        dl = inedges(e,5);
        uL = u(EL,:);
        uR = u(ER,:);  
        
        %compute numerical flux
        [F_hat,Se,error] = HLLE(uL, uR, n, gamma);   

        %increment left residual, decrement right residual
        R(EL,:) = R(EL,:) + F_hat*dl;
        R(ER,:) = R(ER,:) - F_hat*dl;
        
        %add wave speed to tallies
        S(EL) = S(EL) + Se*dl;
        S(ER) = S(ER) + Se*dl;
        
        %calc mach #'s and add to error ind
        ML = mach(uL, gamma);
        MR = mach(uR, gamma);
        e_ind(e,:) = [abs(ML - MR)*dl e];
    end
    
    %% far-field inflow boundaries
    for b = [1 3 4] 
        if error == 1 %if there was an error
            break
        end
        
        edges = bedges{b};
        for e = 1:size(edges,1)
            %extract data: [nA nB nx ny dl Eb]
            Eb = edges(e,6); 
            n = [edges(e,3) edges(e,4)];
            dl = edges(e,5);
            ub = u(Eb,:); 
            
            %compute numerical flux
            [F_hat,Se,error] = HLLE(ub, ufs, n, gamma);   

            %increment boundary element residual
            R(Eb,:) = R(Eb,:) + F_hat*dl;
            
            %add wave speed to tally
            S(Eb) = S(Eb) + Se*dl;
        end
    end
    
    %% Supersonic outflow boundary
    edges = bedges{2};
    for e = 1:size(edges,1)
        if error == 1 %if there was an error
            break
        end
        
        %extract data: [nA nB nx ny dl Eb]
        Eb = edges(e,6); 
        n = [edges(e,3) edges(e,4)];
        dl = edges(e,5);
        ub = u(Eb,:); 

        %compute numerical flux
        [Fbx,Fby,Vel,c,error] = eulerflux(ub,gamma);
        F_hat = Fbx*n(1) + Fby*n(2);
        vel = dot(Vel,n);
        Se = abs(vel)+c; 

        %increment boundary element residual
        R(Eb,:) = R(Eb,:) + F_hat*dl;

        %add wave speed to tally
        S(Eb) = S(Eb) + Se*dl;
    end
    
    %% Solid inviscid wall boundary (Capsule)
    edges = bedges{5};
    Pe = zeros(size(edges,1),1); %Pressure on each capsule edge
    for e = 1:size(edges,1)
        if error == 1 %if there was an error
            break
        end
        
        %extract data: [nA nB nx ny dl Eb]
        nA = edges(e,1);
        nB = edges(e,2);
        Eb = edges(e,6); 
        n = [edges(e,3) edges(e,4)];
        dl = edges(e,5);
        ub = u(Eb,:); 

        %compute numerical flux
        Vel = [ub(2)/ub(1) ub(3)/ub(1)];
        Veltan = Vel - dot(Vel,n).*n;
        magV2 = (norm(Veltan))^2;
        P = (gamma-1)*(ub(4) - 0.5*ub(1)*magV2); %is this right? Yes, want to force the flow to be tangential so calc Pressures with the tangential velocity
        if P < 0 %throw error if pressure is negative 
            error = 1;
            F_hat = 0;
            Se = 0;
        else
            F_hat = [0 P*edges(e,3) P*edges(e,4) 0];
            c = sqrt(gamma*P/ub(1));
            Se = c; 
        end

        %increment boundary element residual
        R(Eb,:) = R(Eb,:) + F_hat*dl;

        %add wave speed to tally
        S(Eb) = S(Eb) + Se*dl;
        
        %calc mach #'s and add to error ind
        Velnorm = dot(Vel,n).*n; %is this right? Velnorm should be zero...? Yes, in the true physicl case the normal velocity would be 0, so the error is calculated with this info.
        unorm = [ub(1) ub(1)*Velnorm(1) ub(1)*Velnorm(2) ub(4)];
        Mb = mach(unorm, gamma);
        Nin = size(inedges,1); %number of interior edges
        e_ind(e+Nin,:) = [abs(Mb)*dl e+Nin];
        
        %% Update Aerodynamic Parameters
        Fx = Fx + P*edges(e,3)*dl; 
        Fy = Fy + P*edges(e,4)*dl;
        
        line = [V(nA,:);V(nB,:)]; %Coordinates of current edge
        r = [mean(line(:,1)) mean(line(:,2))]; %r vector to mid point of current edge
        M = M + (P * (r(1)*n(2)-r(2)*n(1)) * dl);
        
        Pe(e) = P; %Record pressure on each edge
    end
    %% Calculate Aerodynamic Coefficients
    L = Fy*cosd(alphafs) - Fx*sind(alphafs);
    D = Fy*sind(alphafs) + Fx*cosd(alphafs);
    
    %u = [rho, rho*u, rho*v, rho*E];
    Vel_inf = [ufs(2)/ufs(1) ufs(3)/ufs(1)];
    magV2_inf = (norm(Vel_inf))^2;
    rho_inf = ufs(1);
    d = 1.2; %by geometry of capsule
    q = 0.5*rho_inf*magV2_inf; %dymanic pressure of free stream
    P_inf = (gamma-1)*(ufs(4) - 0.5*ufs(1)*magV2_inf);
    
    %Add current coeffs to log
    Cl(N) = L/q*d;
    Cd(N) = D/q*d;
    Cm(N) = M/(q*d^2);
    
    %Overwite current Cp values
    Cp = (Pe - P_inf)./q;
    
    %Reset Fx Fy, and M tallies for new cycle 
    Fx = 0;
    Fy = 0;
    M = 0;
    
    %% Check for convergence
    %MagR_L2 = norm([norm(R(:,1)) norm(R(:,2)) norm(R(:,3)) norm(R(:,4))]); %big 'ol problem with this. Doesnt normalize for number of elements. Going to write by hand.
    SumR = 0;
    nR = size(R,1)*4; %total number of entries in R
    R_list = reshape(R,[nR,1]);
    
    for i = 1:nR
        SumR = SumR + R_list(i)^2; %Cumulative sum of the elements in R_list squared
    end
    
    Rnorm(N) = sqrt(SumR);
    
    if  Rnorm(N) < 10^(-5)
        break
    end
    
    %% Update states
    if error == 1 %if there was an error on this iteration
        CFLTuner = 1; %Back off on the CFL
        u = u_safe; %reset to previous safe state
        continue %Go to next iteration
    end    
    
    u_safe = u; %save state as safe state
    
    CFL = CFLmax - (0.5*CFLTuner); %turn down CFL by 0.5 if CFL tuner is turned on
    for i = 1:size(E,1)
        dt_A = 2*CFL/S(i); %Calculate local timestep/area 
        u(i,:) = u(i,:) - dt_A.*R(i,:);  %update state
    end
    
    %% Check if state output is due
    if count >= count_max
        % Write out restart file
        outputstate = fopen('U.txt','w');
        fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E\n',u');
        fclose(outputstate);
        %reset counter
        count = 0;
        %let user now code is still processing
        fprintf('Working... (N = %d)\n',N);
        
        CFLTuner = 0; %Reset Tuner after count is complete
    end
end
%% Trim Logs
Rnorm = Rnorm(1:N-1);
Cl = Cl(1:N-1);
Cd = Cd(1:N-1);
Cm = Cm(1:N-1);

%% Write out final restart file
outputstate = fopen('U.txt','w');
fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E\n',u');
fclose(outputstate);

end
function [V,E,B,u] = meshadapt(V,E,B,u,e_ind)

% Test Data 1
%{

E = [1,2,4;2,5,4;2,3,5;3,6,5];
V = [0,0;1,0;2,0;0,1;1,1;2,1];
B = cell(4,3);
B{1,1} = 'Left';    B{1,2} = 1; B{1,3} = [4 1];
B{2,1} = 'Top';     B{2,2} = 2; B{2,3} = [6 5; 5 4];
B{3,1} = 'Right';   B{3,2} = 1; B{3,3} = [3 6];
B{4,1} = 'Bottom';  B{4,2} = 2; B{4,3} = [1 2; 2 3];
end
%}

uf = ones(size(E,1),1);

%get boundary egde data
bedges = cell(size(B,1),1);
for j = 1:size(B,1)
    bedges{j} = bedgedat(E,V,B{j,3});
end
%[nA nB nx ny dl Eb flag]

%get interior egde data
C = connect(E);
inedges = inedgedat(E,V,C);
%[nA nB nx ny dl EL ER flag]

%% Flag edges for refinement
Nin = size(inedges,1); %number of interior edges
e_ind = sortrows(e_ind,-1); %sort by decending error 
N_e_ind = size(e_ind,1); %number of entries in e_ind
N_flag = round(0.03*N_e_ind); %small fraction of highest errors to flag
e_ind(N_flag+1:N_e_ind,:) = []; %delete unflagged rows

for e = 1:size(e_ind,1)
    if e_ind(e,2) <= Nin
        inedges(e,8) = 1; %flag interior edge
    elseif e_ind(e,2) > Nin
        bedges{5}(e-Nin,7) = 1; %flag capsule edges
    end    
end
    
%% Set up data structure
    
E_flag = zeros(size(E,1),5); %[is element of row index flagged, number of flagged edges, indicies of edges in e_flag (x3)]
e_flag = zeros(3*size(E,1),3); %point numbers for flagged edges [nA nB nMP] third entry is the midpoint of the first two
flag_count = 0; %counter for number of flagged edges    

%% Mark cell to have all edges flagged if one edge is flagged

%interior edges
for e = 1:size(inedges,1)        
    if inedges(e,8) == 1 %Check if edge is flagged
        flag_count = flag_count + 1; %add to flag count
        
        %Get edge data
        EL = inedges(e,6);
        ER = inedges(e,7);
        
        %Flag elements adjacent to edge, add to their flagged edge tallies,
        %and add index of flagged edge
        E_flag(EL,1) = 1;
        E_flag(ER,1) = 1;  
        E_flag(EL,2) = E_flag(EL,2) + 1; 
        E_flag(ER,2) = E_flag(ER,2) + 1; 
        E_flag(EL,E_flag(EL,2)+2) = flag_count; 
        E_flag(ER,E_flag(ER,2)+2) = flag_count; 
        
        %add new point to V
        nA = inedges(e,1);
        nB = inedges(e,2);

        x = [V(nA,1); V(nB,1)];
        y = [V(nA,2); V(nB,2)];
        mp = [mean(x) mean(y)]; 
        if ~(ismember(mp,V,'rows')) %determine if mp is not already in V
            nV = size(V,1); %number of verticies
            V(nV+1,:) = [mp(1) mp(2)]; %add new vertex to V 
        end
        [~,nMP] = ismember(mp,V,'rows'); %get index of the new mid point in V

        %add nodes to e_flag  
        e_flag(flag_count,1) = nA;
        e_flag(flag_count,2) = nB;
        e_flag(flag_count,3) = nMP;
    end        
end

%boundary edges
for b = 1:size(B,1) 
    edges = bedges{b};
    for e = 1:size(edges,1)
        if edges(e,7) == 1 %Check if edge is flagged
            flag_count = flag_count + 1; %add to flag count
            
            %Get edge data
            Eb = edges(e,6); 

            %Flag elements adjacent to edge, add to their flagged edge tallies,
            %and add index of flagged edge
            E_flag(Eb,1) = 1; 
            E_flag(Eb,2) = E_flag(Eb,2) + 1; 
            E_flag(Eb,E_flag(Eb,2)+2) = flag_count;     

            %add new point to V
            nA = edges(e,1);
            nB = edges(e,2);

            x = [V(nA,1); V(nB,1)];
            y = [V(nA,2); V(nB,2)];
            
            if strcmp(B{b,1},'Capsule') && (x(1) <= 0 || x(2) <= 0)  %special case for heat shield
                disp('Heat Sheild!')
                a1 = atand(y(1)/(.8 - x(1)));
                a2 = atand(y(2)/(.8 - x(2)));
                amp = mean([a1 a2]);
                xmp = -cosd(amp)+0.8;
                ymp = sind(amp);
                mp = [xmp ymp];
                if ~(ismember(mp,V,'rows')) %determine if mp is not already in V
                    nV = size(V,1); %number of verticies
                    V(nV+1,:) = [mp(1) mp(2)]; %add new vertex to V 
                end
                [~,nMP] = ismember(mp,V,'rows'); %get index of the new mid point in V
            else
                mp = [mean(x) mean(y)]; 
                if ~(ismember(mp,V,'rows')) %determine if mp is not already in V
                    nV = size(V,1); %number of verticies
                    V(nV+1,:) = [mp(1) mp(2)]; %add new vertex to V 
                end
                [~,nMP] = ismember(mp,V,'rows'); %get index of the new mid point in V
            end

            %add nodes to e_flag  
            e_flag(flag_count,1) = nA;
            e_flag(flag_count,2) = nB;
            e_flag(flag_count,3) = nMP; 
            
            %add new edges to B and delete original edge
            B_dat = B{b,3};
            nE = size(B_dat,1); %nuber of edges in B_dat
            B_dat(nE+1,:) = [nB nMP];
            B_dat(nE+2,:) = [nMP nA];
            
            [~,index_del] = ismember([nB nA],B_dat,'rows');
            B_dat(index_del,:) = [];
            B{b,3} = B_dat;
            
        end 
    end
end  

%% Flag all other edges for flagged cells

for index_fE = find(E_flag(:,1))' %find element indicies that are flagged
    
    %interior edges 
    e = find((inedges(:,6) == index_fE) | (inedges(:,7) == index_fE)); %find all interior edges that boarder a flagged cell
    for idx_e = e'
        if inedges(idx_e,8) == 0 %if the edge hasnt been flagged, proceede  
            inedges(idx_e,8) = 1; %flag the edge   
            flag_count = flag_count + 1; %add to flag count 

            %Flag elements adjacent to edge, add to their flagged edge tallies,
            %and add index of flagged edge
            EL = inedges(idx_e,6); 
            ER = inedges(idx_e,7);     
            E_flag(EL,1) = 1;
            E_flag(ER,1) = 1;  
            E_flag(EL,2) = E_flag(EL,2) + 1;
            E_flag(ER,2) = E_flag(ER,2) + 1; 
            E_flag(EL,E_flag(EL,2)+2) = flag_count;
            E_flag(ER,E_flag(ER,2)+2) = flag_count;
            
            %add new point to V
            nA = inedges(idx_e,1);
            nB = inedges(idx_e,2);
            
            x = [V(nA,1); V(nB,1)];
            y = [V(nA,2); V(nB,2)];
            mp = [mean(x) mean(y)]; 
            if ~(ismember(mp,V,'rows')) %determine if mp is not already in V
                nV = size(V,1); %number of verticies
                V(nV+1,:) = [mp(1) mp(2)]; %add new vertex to V 
            end
            [~,nMP] = ismember(mp,V,'rows'); %get index of the new mid point in V
            
            %add nodes to e_flag  
            e_flag(flag_count,1) = nA;
            e_flag(flag_count,2) = nB;
            e_flag(flag_count,3) = nMP;
        end
    end
    
    %boundary edges
    for b = 1:size(B,1)
        edges = bedges{b};
        e = find(edges(:,6) == index_fE); %find all boundary edges that boarder a flagged cell
        for idx_e = e' 
            if edges(idx_e,7) == 0 %if the edge hasnt been flagged, proceede 
                edges(idx_e,7) = 1; %flag the edge
                flag_count = flag_count + 1; %add to flag count 

                %Flag elements adjacent to edge, add to their flagged edge tallies,
                %and add index of flagged edge
                Eb = edges(e,6); 
                E_flag(Eb,1) = 1; 
                E_flag(Eb,2) = E_flag(Eb,2) + 1; 
                E_flag(Eb,E_flag(Eb,2)+2) = flag_count;
                
                
                %add new point to V
                nA = edges(idx_e,1);
                nB = edges(idx_e,2);

                x = [V(nA,1); V(nB,1)];
                y = [V(nA,2); V(nB,2)];
                
                if strcmp(B{b,1},'Capsule') && (x(1) <= 0 || x(2) <= 0)  %special case for heat shield
                disp('Heat Sheild!')
                a1 = atand(y(1)/(.8 - x(1)));
                a2 = atand(y(2)/(.8 - x(2)));
                amp = mean([a1 a2]);
                xmp = -cosd(amp)+0.8;
                ymp = sind(amp);
                mp = [xmp ymp];
                if ~(ismember(mp,V,'rows')) %determine if mp is not already in V
                    nV = size(V,1); %number of verticies
                    V(nV+1,:) = [mp(1) mp(2)]; %add new vertex to V 
                end
                [~,nMP] = ismember(mp,V,'rows'); %get index of the new mid point in V
                else
                    mp = [mean(x) mean(y)]; 
                    if ~(ismember(mp,V,'rows')) %determine if mp is not already in V
                        nV = size(V,1); %number of verticies
                        V(nV+1,:) = [mp(1) mp(2)]; %add new vertex to V 
                    end
                    [~,nMP] = ismember(mp,V,'rows'); %get index of the new mid point in V
                end

                %add nodes to e_flag  
                e_flag(flag_count,1) = nA;
                e_flag(flag_count,2) = nB;
                e_flag(flag_count,3) = nMP;  

                %add new edges to B and delete original edge
                B_dat = B{b,3};
                nE = size(B_dat,1);
                B_dat(nE+1,:) = [nB nMP];
                B_dat(nE+2,:) = [nMP nA];

                [~,index_del] = ismember([nB nA],B_dat,'rows');
                B_dat(index_del,:) = [];
                B{b,3} = B_dat;
            end
        end
    end
end

e_flag = e_flag(1:flag_count,:); %trim e_flag with flag_count

%% Plot Pre-Split
%{
uf2 = uf + E_flag; %mark flagged cells for plotting

%Visualize Flagged Edges and Cells
f = figure(10);
clf(10)
f.OuterPosition = [100 278 560*1.25 420];

subplot('Position',[.05 .125 0.45 .825]);
meshplot(uf2,E,V,'Flag','Flagged Mesh',1,0)
colorbar('off')
hold on

%plot flagged edges
for i = 1:size(e_flag,1)
    x = [V(e_flag(i,1),1); V(e_flag(i,2),1)];
    y = [V(e_flag(i,1),2); V(e_flag(i,2),2)];
    plot(x,y,'r','LineWidth',2);
end
%}

%% Split Flagged Elements

for index_fE = find(E_flag(:,1))' %find element indicies that are flagged
    
    %% one edge flagged
    if E_flag(index_fE,2) == 1
        
        edge = e_flag(E_flag(index_fE,3),:); %get info of flagged edge

        %Get points of original flagged element
        Pts = E(index_fE,:);
        
        %Determine points for upper split
        p1 = edge(3); %new point
        p2 = Pts((Pts ~= edge(1)) & (Pts ~= edge(2))); %point of original element not on flagged edge
        
        %pick other points in original element
        for i = Pts((Pts == edge(1)) | (Pts == edge(2)))          
            if isareapos(V,p1,p2,i)    
                p3 = i; %point of original element that makes area positive
            elseif ~isareapos(V,p1,p2,i)
                pL2 = i; %point on edge for lower split
            end           
        end
        
        nE = size(E,1); %number of elements 
         
        %Add upper split element with corresponding state
        E(nE+1,:) = [p1 p2 p3];
        u(nE+1,:) = u(index_fE,:);
        
        %Add lower split element with corresponding state
        p1 = edge(3); %new point
        p3 = p2; %point of original element not on flagged edge
        p2 = pL2; 
        E(nE+2,:) = [p1 p2 p3];
        u(nE+2,:) = u(index_fE,:);
        
     end
    
    %% two edges flagged
    if E_flag(index_fE,2) == 2
        
        %get info of flagged edges and add new points to V
        edge = zeros(2,3);
        for e = 1:size(edge,1)
            edge(e,:) = e_flag(E_flag(index_fE,e+2),:);
        end
        
        %Get points of original flagged element
        Pts = E(index_fE,:);
        
        nE = size(E,1); %number of elements
         
        %determine points of element ajacent to both flagged edges
        p1 = intersect(edge(1,:),edge(2,:)); %start at shared point
        
        %find combo of other two mp on edges 1 and 2 that makes area positive
        if isareapos(V,p1,edge(1,3),edge(2,3))
            p2 = edge(1,3);
            p3 = edge(2,3);
        else
            p3 = edge(1,3);
            p2 = edge(2,3);
        end
        
        %Add double adjacent element with corresponding state
        E(nE+1,:) = [p1 p2 p3];
        u(nE+1,:) = u(index_fE,:);
        
        %Calculate the angles between each flagged edge and the nonflagged edge    
        edge_nf = setdiff(Pts,p1); %non flagged edge
        
        r1 =  V(edge(1,2),:) - V(edge(1,1),:); %vector for edge 1
        r2 =  V(edge(2,2),:) - V(edge(2,1),:); %vector for edge 2
        r_nf = V(edge_nf(2),:) - V(edge_nf(1),:); %vector for non flagged edge
        
        a1 = acosd(dot(r1,r_nf)/(norm(r1)*norm(r_nf)));
        a2 = acosd(dot(r2,r_nf)/(norm(r2)*norm(r_nf)));
        
        %slpit along larger angle
        if a1 >= a2
            p1 = intersect(edge(1,1:2),edge_nf);
            
            %find combo of other two mp that makes area positive
            if isareapos(V,p1,edge(1,3),edge(2,3))
                p2 = edge(1,3);
                p3 = edge(2,3);
            else
                p3 = edge(1,3);
                p2 = edge(2,3);
            end
            
            %add edge 1 adjacent element with corresponding state
            E(nE+2,:) = [p1 p2 p3];
            u(nE+2,:) = u(index_fE,:);
            
            %find combo of edge 2 mp and the other edge_nf point that makes area positive
            if isareapos(V,p1,intersect(edge(2,1:2),edge_nf),edge(2,3))
                p2 = intersect(edge(2,1:2),edge_nf);
                p3 = edge(2,3);
            else
                p3 = intersect(edge(2,1:2),edge_nf);
                p2 = edge(2,3);
            end
            
            %add edge 2 adjacent element with corresponding state
            E(nE+3,:) = [p1 p2 p3];
            u(nE+3,:) = u(index_fE,:);
        elseif a1 < a2
            p1 = intersect(edge(2,1:2),edge_nf);
            
            %find combo of other two mp that makes area positive
            if isareapos(V,p1,edge(1,3),edge(2,3))
                p2 = edge(1,3);
                p3 = edge(2,3);
            else
                p3 = edge(1,3);
                p2 = edge(2,3);
            end
            
            %add edge 2 adjacent element with corresponding state
            E(nE+2,:) = [p1 p2 p3];
            u(nE+2,:) = u(index_fE,:);
            
            %find combo of edge 1 mp and the other edge_nf point that makes area positive
            if isareapos(V,p1,intersect(edge(1,1:2),edge_nf),edge(1,3))
                p2 = intersect(edge(1,1:2),edge_nf);
                p3 = edge(1,3);
            else
                p3 = intersect(edge(1,1:2),edge_nf);
                p2 = edge(1,3);
            end
            
            %add edge 2 adjacent element with corresponding state
            E(nE+3,:) = [p1 p2 p3];
            u(nE+3,:) = u(index_fE,:);
        end      
    end
    
    %% three edges flagged
    if E_flag(index_fE,2) == 3
        
        %get info of flagged edges and add new points to V
        edge = zeros(3,3);
        for e = 1:size(edge,1)
            edge(e,:) = e_flag(E_flag(index_fE,e+2),:);
        end
      
        %Determine points for central split
        p1 = edge(1,3); %pick a starting mp
        
        %find combo of other two mp that makes area positive
        if isareapos(V,p1,edge(2,3),edge(3,3))
            p2 = edge(2,3);
            p3 = edge(3,3);
        else
            p3 = edge(2,3);
            p2 = edge(3,3);
        end
        
        nE = size(E,1); %number of elements
         
        %Add central split element with corresponding state
        E(nE+1,:) = [p1 p2 p3];
        u(nE+1,:) = u(index_fE,:);
        
        %determine points for element along edge 1 and edge 2 (element 1,2)
        p1 = intersect(edge(1,:),edge(2,:)); %start at shared point
        
        %find combo of other two mp on edges 1 and 2 that makes area positive
        if isareapos(V,p1,edge(1,3),edge(2,3))
            p2 = edge(1,3);
            p3 = edge(2,3);
        else
            p3 = edge(1,3);
            p2 = edge(2,3);
        end
        
        %Add element 1,2 with corresponding state
        E(nE+2,:) = [p1 p2 p3];
        u(nE+2,:) = u(index_fE,:);
        
        %determine points for element along edge 2 and edge 3 (element 2,3)
        p1 = intersect(edge(2,:),edge(3,:)); %start at shared point
        
        %find combo of other two mp on edges 2 and 3 that makes area positive
        if isareapos(V,p1,edge(2,3),edge(3,3))
            p2 = edge(2,3);
            p3 = edge(3,3);
        else
            p3 = edge(2,3);
            p2 = edge(3,3);
        end
        
        %Add element 2,3 with corresponding state
        E(nE+3,:) = [p1 p2 p3];
        u(nE+3,:) = u(index_fE,:);
        
        %determine points for element along edge 1 and edge 3 (element 1,3)
        p1 = intersect(edge(1,:),edge(3,:)); %start at shared point
        
        %find combo of other two mp on edges 1 and 3 that makes area positive      
        if isareapos(V,p1,edge(1,3),edge(3,3))
            p2 = edge(1,3);
            p3 = edge(3,3);
        else
            p3 = edge(1,3);
            p2 = edge(3,3);
        end
        
        %Add element 1,3 with corresponding state
        E(nE+4,:) = [p1 p2 p3];
        u(nE+4,:) = u(index_fE,:);
        
    end 
end

%delete original flagged elements and corresponding state
index_fE = find(E_flag(:,1));
E(index_fE,:) = []; 
u(index_fE,:) = [];

%% Plot after split

%{
subplot('Position',[.525 .125 0.45 .825]); %[left bottom width height]

%Visualize Flagged Edges and Cells
uf = ones(size(E,1),1);
meshplot(uf,E,V,'Ones','Refined Mesh',1,0)
colorbar('off')
hold on
%scatter(V(:,1),V(:,2),'k','filled') 

%get boundary egde data
bedges = cell(size(B,1),1);
for j = 1:size(B,1)
    bedges{j} = bedgedat(E,V,B{j,3});
end

%plot boundary edges
for j = 1:size(bedges,1)
    edge = bedges{j};
    c = 0;
    for i = 1:size(edge,1)
        if c == 1
            x = [V(edge(i,1),1); V(edge(i,2),1)];
            y = [V(edge(i,1),2); V(edge(i,2),2)];
            plot(x,y,'c','LineWidth',2);
            c = c*-1;
        elseif c == -1
            x = [V(edge(i,1),1); V(edge(i,2),1)];
            y = [V(edge(i,1),2); V(edge(i,2),2)];
            plot(x,y,'r','LineWidth',2);
            c = c+1;
        elseif c == 0
            x = [V(edge(i,1),1); V(edge(i,2),1)];
            y = [V(edge(i,1),2); V(edge(i,2),2)];
            plot(x,y,'y','LineWidth',2);
            c = c+1;
        end
    end
end

%get interior egde data
C = connect(E);
inedges = inedgedat(E,V,C);

%plot interior edges
for i = 1:size(inedges,1)
    x = [V(inedges(i,1),1); V(inedges(i,2),1)];
    y = [V(inedges(i,1),2); V(inedges(i,2),2)];
    plot(x,y,'k','LineWidth',0.1);
end

input('Done...?');
fprintf('--------> Adaptation Plots: Done!\n\n');
%}

%% Write New Mesh File
fname = input('Enter new mesh file save name: ','s');
fileID = fopen(fname,'w');
nNode = size(V,1);
nElem = size(E,1);
Dim = 2;
Nb = size(B,1);

fprintf(fileID,'%d %d %d\n',nNode,nElem,Dim);
fprintf(fileID,'%E %E\n',V');

fprintf(fileID,'%d\n',Nb);

for i = 1:Nb
    nBface = size(B{i,3},1);
    Title = char(B{i,1});
    fprintf(fileID,'%d %d %s\n',nBface,2,Title);
    fprintf(fileID,'%d %d\n',B{i,3}');
end  
fprintf(fileID,'%d %d %s\n',nElem,1,'TriLagrange');
fprintf(fileID,'%d %d %d\n',E');

fclose(fileID);

%% Write out restart file
expU = input('Enter new expanded state file name: ','s');
outputstate = fopen(expU,'w');
fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E\n',u');
fclose(outputstate);

end
    
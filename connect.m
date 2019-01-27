function C = connect(E)
% C = connect(E) determines the shared edges of the elements E of the mesh
% 
%   C is of the form C = [E1 e1 E2 e2]; where E2 and E2 refer to the two
%   elements that share a connection, and e1 and e2 are the respective
%   local edge numbers where edge 1 is across from point 1 in E.

%Determine size of sparse matrix and create
sizeS = max(max(E));
S = spalloc(sizeS,sizeS,6*length(E(:,1))); %6*length defines amount of non zero elements to be stored. guess average of 6 connections per element 

%Create output matrix 
C = zeros(6*length(E(:,1)),4);
c = 1; %arbitrary counter used to trim C at end of loop

%Loop over all elements
for i = (1:length(E(:,1))) %i = triange number
    
    %loop over edges
    for j = (1:3) %all triangles so always 3 edges; j = edge number
        
        e = E(i,:); %restore triangle i
        e(j) = []; %eliminate point j to define edge as array of other two points
        e = sort(e); %sort points to correctly ID location in S
        
        if full(S(e(1),e(2))) ~= 0 %exists entry at edge location
           
            t0 = full(S(e(1),e(2))); %original triangle number 
            
            %need original edge number
            %E(t1,:) is the original element
            %e(1) and e(2) are sequential point numbers that exist in E(t1)
            
            e0 = find(E(t0,:) ~= e(1) & E(t0,:) ~= e(2)); %index of point that is not e(1) or e(2)
            
            if t0 < i %correctly order t1 and t2
                t1 = t0;
                e1 = e0;
                t2 = i;
                e2 = j;
            else
                t1 = i;
                e1 = j;
                t2 = t0;
                e2 = e0;
            end
            
            C(c,:) = [t1 e1 t2 e2]; %create new row in C
            c = c+1; %increment c
            
        else
            %store triange and edge number in S
            S = S + sparse(e(1),e(2),i,sizeS,sizeS,6*length(E(:,1))); 
            %S(e(1),e(2)) = i
            
        end
    end
end

C = C(1:c-1,:); %trim C with c

C = sortrows(C,[1 3]); %sort matrix C with correct rules

end

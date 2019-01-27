function inedges = inedgedat(E,V,C)
% inedgedat gives the interior edge data for the mesh
%   
%   the output (inedges) is a matrix of the form [nA nB nx ny dl EL ER flag] with 
%   length equal the number of interior edges, where:
%   nA is the index in V of node A
%   nB is the index in V of node B
%   nx and ny are the components of the unit normal vector respectively
%   (unit normal points from L to R)
%   dl is the length of the edge
%   EL is the index of the left element in E
%   ER is the index of the right element in E
%   flag is either 1 or 0 to indicate edge for refinement

inedges = zeros(size(C,1),8);

for i = 1:size(C,1)
    EL = C(i,1); ER = C(i,3);
    tri = E(EL,:); %gives three node indicies of EL
    
    %conditions ensure that nA is always on top, and therefor normal points
    %towards R
    if C(i,2) == 1
        nA = tri(3); nB = tri(2);
    elseif C(i,2) == 2
        nA = tri(1); nB = tri(3);
    elseif C(i,2) == 3
        nA = tri(2); nB = tri(1);
    end
    
    VA = V(nA,:); VB = V(nB,:); %coordinates for both nodes
    
    dl = ((VA(1)-VB(1))^2 + (VA(2)-VB(2))^2)^0.5;
    
    nx = (VA(2)-VB(2))/dl; ny = (VB(1)-VA(1))/dl;
    
    %assign values to matrix
    inedges(i,1) = nA;
    inedges(i,2) = nB;
    inedges(i,3) = nx;
    inedges(i,4) = ny;
    inedges(i,5) = dl;
    inedges(i,6) = EL;
    inedges(i,7) = ER;
end
end


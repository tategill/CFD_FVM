function bedges = bedgedat(E,V,B)
% dedgedat gives the boundary edge data for the boundary group given the list of
% nodes in the boundary group. (nodes in B should be listed top to bottom)
%   
%   the output (bedges) is a matrix of the form [nA nB nx ny dl Eb flag] with 
%   length equal the number of boundary edges in the group, where:
%   nA is the index in V of node A
%   nB is the index in V of node B
%   nx and ny are the components of the unit normal vector respectively
%   (unit normal points outward from numerical domain L)
%   dl is the length of the edge
%   Eb is the index of the boundary element in E
%   flag is either 1 or 0 to indicate edge for refinement

bedges = zeros(size(B,1),7);

for i = 1:size(B,1)
    nA = B(i,2); nB = B(i,1);
    
    Eb = find(and(any(E == nA, 2), any(E == nB, 2)));   
    
    VA = V(nA,:); VB = V(nB,:); %coordinates for both nodes
    
    dl = ((VA(1)-VB(1))^2 + (VA(2)-VB(2))^2)^0.5;
    
    nx = (VA(2)-VB(2))/dl; ny = (VB(1)-VA(1))/dl;
    
    %assign values to matrix
    bedges(i,1) = nA;
    bedges(i,2) = nB;
    bedges(i,3) = nx;
    bedges(i,4) = ny;
    bedges(i,5) = dl;
    bedges(i,6) = Eb;
  
end
end


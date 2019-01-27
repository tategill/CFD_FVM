function [V, E, B] = readgri(fname)
% readgri returns the mesh data in the fname 
%
%   V is the Vertex list: V = [x1 y1; x2 y2; ect..]
%
%   E is the Elememt list: E = [p1 p2 p3; p1 p2 p3; ect..] where p's refer to
%   points via row numbers in V
%
%   B is a cell array that contains the boundary titles in the first column
%   the number of edges in the second column and the points refering to
%   those edges in the third column (in a similar fashion to E).

    f = fopen(fname, 'r');
    
    %read and assign mesh info
    infoln = str2double(split(fgetl(f)));
    Nn = infoln(1); Ne = infoln(2); Dim = infoln(3); 
    
    %read vertices
    V = fscanf(f,'%f',[Dim,Nn])';
    
    %read boundaries
    Nb = fscanf(f,'%f',[1,1]);
    B = cell(Nb,3);
    for i = 1:Nb
        Bgroupi = textscan(f,'%f %f %s',[1,3]);
        B{i,1} = Bgroupi{3}; 
        B{i,2} = Bgroupi{1}; 
        B{i,3} = fscanf(f,'%f',[Bgroupi{2},Bgroupi{1}])';
    end  
    
    %read elements
    Ne0 = 0; E = zeros(Ne,3);
    while Ne0 < Ne 
        Einfo = textscan(f,'%f %f %s',[1,3]);
        ne = Einfo{1};
        Ei = fscanf(f,'%f',[3, ne])';
        E(Ne0+1:ne,:) = Ei;
        Ne0 = Ne0 + ne;   
    end
    
    fclose(f);
    
end
    

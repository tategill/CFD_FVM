function plotcapsule(B,V)

capB = sortrows(B{5,3},1);
x = zeros(size(capB,1),1);
y = x;

for i = 1:size(capB,1)
    x(i) = V(capB(i,1),1);
    y(i) = V(capB(i,1),2);
end

pts = sortrows([x y],[1 2]);
ptsUp = pts(pts(:,2)>0,:);
ptsUp = sortrows(ptsUp,[1 -2]);
ptsDn = flip(pts(pts(:,2)<0,:));

%{
% for i = 1:size(ptsUp,1)
%     %input('next pt')
%     scatter(ptsUp(i,1),ptsUp(i,2),'k');
% end

% for i = 1:size(ptsDn,1)
%     %input('next pt')
%     scatter(ptsDn(i,1),ptsDn(i,2),'k');
% end
%}

x = [ptsUp(:,1);ptsDn(:,1)];
y = [ptsUp(:,2);ptsDn(:,2)];

fill(x,y,[0.35 0.35 0.35])



end
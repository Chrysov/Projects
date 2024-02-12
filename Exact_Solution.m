function E=Exact_Solution(k)
%We will follow the assignment hints to produce the exact mesh.
%Introduction of h, and x and y ranges.

h=1/(k+1);
y=(0:h:1);
x=(0:h:1);
[X,Y]=meshgrid(x,y);

%Set exact solution which is given.
u=sin(2*pi*X).*sin(2*pi*Y)+10*(X.^4-X.^2).*(Y.^4-Y.^2)+2;

%Mesh of the exact solution.
r=mesh(X,Y,u);
title('U^h');
r.FaceColor='Flat';
end

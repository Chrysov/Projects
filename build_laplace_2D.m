function A=build_laplace_2D(k)

%Arrays of 1's that we will use.
e1=ones(k^2,1);
e2=ones(k,1);

%Creating our block matrix T.
D1= spdiags([e2,-4*e2,e2], -1:1, k, k);
%Making a repmat k^2 by k^2 matrix and keeping only its three central diagonals.
J=repmat(D1,k,k);
full(J)
x1=diag(diag(J,-1),-1);
x2=diag(diag(J));
x3=diag(diag(J,1),1);
full(x1)
full(x2)
full(x3)
X=x1+x2+x3;

%Our I matrices on our off diagonals take the form.
D2=spdiags([e1 e1], [-k k], k^2, k^2);

%Adding the two matrices and multiplying them by 1/(h^2)
A=D2+X;
h = 1/(k+1);
A = A/h^2;
spy(A)
end

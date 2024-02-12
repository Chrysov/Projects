function P=Parabolic(k)

%Given values
eta=2;
alpha=5;
to=0.05;

%Spacial Info
dx=30/(k+1);
x=(-10:dx:20);

%Given function u at time to
u=(1/(to^(1/2)))*exp((-((x-alpha*to).^2))/(4*eta*to));
u=reshape(u,[k+2,1]);

%2 Different Timesteps.
Dtd= (dx^2)/(2*eta);
Dta=(dx)/(alpha);

%(1) EXPLICIT METHOD

%Given timestep
Dt1=(0.4750)*min(Dta,Dtd);

%New timestep
Dtnew=((1-to)/(ceil((1-to)/(Dt1))));

%Now we create our matrix with spdiags like in Q6.

e=ones(k+2,1);
Matr1= spdiags([e*(((eta*Dtnew)/(dx^2))+((alpha*Dtnew)/(dx))),...;
    e*(1-((2*eta*Dtnew)/(dx^2))-((alpha*Dtnew)/(dx))),...;
    e*((eta*Dtnew)/(dx^2))], -1:1, k+2, k+2);

%Preparing for our loop
T=1;
tI=0.05;
UL=u;

%Loop
while tI<1
    T=T+1;
    tI=tI+Dtnew;
    Ur=UL(:,end);
    UL=Matr1*Ur;
end

%Exact Plot
figure
EXACT=exp(-(((x-alpha).^2)/(4*eta)));
plot(x,EXACT);
hold on

%Explicit Plot at t=3
plot(x,UL)
hold on

%(2) SEMI-IMPLICIT METHOD

%Given and New Timesteps
Dt2=0.95*Dta;
Dt2new=(1-to)/(ceil((1-to)/(Dt2)));

%Now we have two matrices for which we use spdiags again.
Matr21=spdiags([-e*((eta*Dt2)/(dx^2)),e*(1+((2*eta*Dt2)/(dx^2))),...;
    e*((-eta*Dt2)/(dx^2))],-1:1, k+2, k+2);
Matr22=spdiags([e*((alpha*Dt2)/(dx)),e*(1-((alpha*Dt2)/(dx)))], -1:0, k+2, k+2);

%Preparing for our loop
UR=u;
TR=0.05;

%Our loop
while TR<1
   TR= TR+Dt2;
   UR=Matr22*UR;
   UR=Matr21\UR;   
end

%Implicit Plot
plot(x,UR)
legend('Exact','Explicit','Semi-Implict')
hold off
end










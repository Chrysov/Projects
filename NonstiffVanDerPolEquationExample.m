%% Solve Nonstiff Equation
% The van der Pol equation is a second order ODE
%
% $$y''_1 - \mu \left( 1 - y_1^2\right) y'_1+y_1=0,$$
%
% where $\mu > 0$ is a scalar parameter. Rewrite this equation as a system
% of first-order ODEs by making the substitution $y'_1 = y_2$. The
% resulting system of first-order ODEs is
%
% $$
% \begin{array}{cl}
% y'_1 &= y_2\\
% y'_2 &= \mu (1-y_1^2) y_2 - y_1.\end{array}
% $$
%

%% 
% The function file |vdp1.m| represents the van der Pol equation using $\mu
% = 1$. The variables $y_1$ and $y_2$ are the entries |y(1)| and |y(2)| of
% a two-element vector, |dydt|.
%
% <include>vdp1.m</include>
%

%% 
% Solve the ODE using the |ode45| function on the time interval |[0 20]|
% with initial values |[2 0]|. The resulting output is a column vector of
% time points |t| and a solution array |y|. Each row in |y| corresponds to
% a time returned in the corresponding row of |t|. The first column of |y|
% corresponds to $y_1$, and the second column to $y_2$.
[t,y] = ode45(@vdp1,[0 20],[2; 0]);

%%
% Plot the solutions for $y_1$ and $y_2$ against |t|.
plot(t,y(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')


%% 
% Copyright 2012 The MathWorks, Inc.
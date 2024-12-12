A = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];

b = [2 1 2 1 -1 2];

Aeq = [1 1/4];
beq = 1/2;

lb = [-1,-0.5];
ub = [1.5,1.25];

% cost/objective fn
f = [-1 -1/3];

% switch to 'interior-point'
options = optimoptions('linprog','Algorithm','interior-point');

% solves `x`; min. f'*x and Ax <= B
x = linprog(f,A,b,Aeq,beq,lb,ub,options)


%%
%BONUS
%optimization in each step

clear all;
clc;

N=1000;
h=0.01;
v=1;
wmin=-pi*4;
wmax=pi*4;

%initial conditions
y0=1.5;
h0=0;
r=1;
Q=[1 0; 0 1];
%%
%x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0)
%objective function
H = sparse(eye(N)*r); %identity matrix with size N with diagonal r
for i = 1:N
 H = blkdiag(H, Q);
end
H = 2 * H; 
f = zeros(3*N,1);
%%
%inequality constraint
%Aineq: Matrix in linear inequality constraints Aineq*x ? bineq
%bineq: Vector in linear inequality constraints Aineq*x ? bineq
Aineq = sparse([eye(N) zeros(N,2*N);-1*eye(N) zeros(N,2*N)]);
bineq = sparse([ones(N,1)*(-wmin); ones(N,1)*wmax]);
%%
%equality constraint
%Aeq: Matrix in linear equality constraints Aeq*x = beq
%beq: Vector in linear equality constraints Aeq*x = beq


LEFTMATRIX = zeros(2*N,N);
for i = 1:N
    LEFTMATRIX(2*i,i)=-h;
end

RIGHTMATRIX = zeros(2*N,2*N);
for i = 1:2*N
    if mod(i,2) == 1 %if odd
        RIGHTMATRIX(i,i)=1;
            if (i<2*N-1)%to prevent out of indexing
            RIGHTMATRIX(i+2,i)=-1;
            end
    else %if even
        if i<2*N-1
            RIGHTMATRIX(i,i)=1;
            RIGHTMATRIX(i+1,i)=-h*v;
            RIGHTMATRIX(i+2,i)=-1;
        elseif i<2*N
            RIGHTMATRIX(i,i)=1;
            RIGHTMATRIX(i+1,i)=-h*v;
        else %i=2*N
            RIGHTMATRIX(i,i)=1;
        end
    end
end


Aeq = [LEFTMATRIX RIGHTMATRIX];
beq = [y0+h*v*h0; h0; zeros(2*N-2,1)];
%%
%allocating memory for vectors
y=zeros(1,N+1);
head=zeros(1,N+1);
y(1)=y0;
head(1)=h0;
w=zeros(1,N+1);

for i=2:N+1 % must start at 2 since x(i-1)=x(1) for first iteration
    
    %update initial conditions
    y0=y(i-1);
    h0=head(i-1);
    
    %update constraints that depend on initial conditions
    beq = [y0+h*v*h0; h0; zeros(2*N-2,1)];
    
    %optimize
    [x, fval] = quadprog(H,f,Aineq,bineq,Aeq,beq);
    
    %input w0 value
    w(i-1)=x(1);
 
    %dynamics...y index starts at 1 (initial condition)
    y(i)=y(i-1)+h*v*head(i-1);
    head(i)=head(i-1)+h*w(i-1)+h*pi/25;
end

t=(0:N)*h;
subplot(3,1,1);
plot(t,y);
title('{yk} with Model Predictive Control');
xlabel('Time (s)')
ylabel('Position');

subplot(3,1,2);
plot(t, head);
title('{hk} with Model Predictive Control');
xlabel('Time (s)');
ylabel('Heading');

subplot(3,1,3);
plot(t, w);
title('Turn Rate with Model Predictive Control');
xlabel('Time (s)');
ylabel('w (rad/s)');


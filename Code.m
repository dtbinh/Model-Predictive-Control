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

%x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0)
%objective function
H = sparse(eye(N)*r); %identity matrix with size N with diagonal r
for i = 1:N
 H = blkdiag(H, Q);
end
H = 2 * H; 
f = zeros(3*N,1);

%inequality constraint
%Aineq: Matrix in linear inequality constraints Aineq*x ? bineq
%bineq: Vector in linear inequality constraints Aineq*x ? bineq
Aineq = sparse([eye(N) zeros(N,2*N);-1*eye(N) zeros(N,2*N)]);
bineq = sparse([ones(N,1)*(-wmin); ones(N,1)*wmax]);

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

[x, fval] = quadprog(H,f,Aineq,bineq,Aeq,beq);
% x = [w0 w1 ... wn x1 x2 ...xn]

%construct {w}
w = x(1:N);

yk = zeros(1,1000); %initializing
hk = zeros(1,1000); 

%construct {yk}
loopcounter=1;
for i=N+1:3*N
    if mod(i,2)==1
        yk(loopcounter)= x(i);
        loopcounter = loopcounter +1;
    end
end

%construct {hk}

loopcounter1=1;
for i=N+1:3*N
    if mod(i,2)==0
        hk(loopcounter1)= x(i);
        loopcounter1 = loopcounter1 +1;
    end
end
t=(0:N)*h;
%subplot(3,1,1);
%plot(t,[y0 yk]);
%title('{yk}');
%xlabel('Time (s)')
%ylabel('Position');

%subplot(3,1,2);
%plot(t, [h0, hk]);
%title('{hk}');
%xlabel('Time (s)');
%ylabel('Heading');

%subplot(3,1,3);
%plot(t(1:N), w);
%title('Turn Rate');
%xlabel('Time (s)');
%ylabel('w (rad/s)');
%%
%f) simulation with noise


y=zeros(1,N+1);
head=zeros(1,N+1);
y(1)=1.5;
head(1)=0;

for i=1:N
   y(i+1)=y(i)+h*v*head(i);
   head(i+1)=head(i)+h*w(i)+h*pi/25;
end

t=(0:N)*h;
%plot(t,y); 
%hold on;
%plot(t,[y0 yk]);
%legend('with noise', 'without noise');
%xlabel('time (s)');
%ylabel('Position');
%title('{yk}');
plot(t,head);
hold on;
plot(t,[h0 hk]);
legend('with noise', 'without noise');
xlabel('time (s)');
ylabel('Heading');
title('{hk}');











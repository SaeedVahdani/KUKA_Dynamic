clear all
close all
clc

syms q1 q2 q3 q4 real
syms q_d1 q_d2 q_d3 q_d4 real
Links = 4;
Type = 'PRRR+';
g = [0,-9.81,0];
q = [q1, q2, q3, q4]';
q_d = [q_d1, q_d2, q_d3, q_d4]';

%% Base Properties
% mass{1} = 40.117;
% CoM{1} = [0,-0.033,0];
% MoI{1} = [13.366,13.428,0.110]; 
% PoI{1} = [0,0,0];              
R{1} = eye(3);
P{1} = [0,0,0];

%% Link 1 Properties
mass{1} = 22.064;
CoM{1} = [-0.031,0.094,0];
MoI{1} = [0.134,0.165,0.186]; 
PoI{1} = [0,0.001,-0.012];        
R{2} = [1 0 0; 0 0 1;0 -1 0];
P{2} = [0,0.197,0];

%% Link 2 Properties
mass{2}= 9.805;
CoM{2}= [-0.005,-0.001,0.082];
MoI{2}= [0.082,0.053,0.058]; 
PoI{2}= [-0.001,0,0];        
R{3}= [1 0 0; 0 0 -1; 0 1 0]*[0 -1 0;1 0 0; 0 0 1];
P{3}= [-0.00059,0,0.1522];

%% Link 3 Properties
mass{3}= 13.302;
CoM{3}= [0.118,0.011,-0.001];
MoI{3}= [0.048,0.175,0.153]; 
PoI{3}= [0,0.001,0.003];        
R{4}= [0 1 0; -1 0 0; 0 0 1];
P{4}= [0.290,0,0];

%% Link 4 Properties
mass{4}= 10.723;
CoM{4}= [0.143,0.015,-0.001];
MoI{4}= [0.018,0.166,0.164]; 
PoI{4}= [0,-0.001,0.009];        
R{5}= eye(3);
P{5}= [0.3826,-0.020,0];

%% Transformation Matrix
A = sym(zeros(4,4,5));


for i = 1:5
    A(:,:,i) = [R{i},   P{i}'  ; 0,0,0,1];
    
    if Type(i)=='P'
       A(:,:,i) = A(:,:,i)*[eye(3), [0;0;q(i)] ; 0,0,0,1];
    elseif Type(i)=='R'
            R2= [cos(q(i)),    -sin(q(i))    ,0
                 sin(q(i)),    cos(q(i))     ,0
                 0     ,    0          ,1];
            T = [R2, [0;0;0]  ; 0,0,0,1];
            A(:,:,i) = A(:,:,i)*T;
    end
end


T(:,:,1) = A(:,:,1);
for i = 2:5
   T(:,:,i) =  T(:,:,i-1)*A(:,:,i);
end


%% Inertia Matrix
I = zeros(4,4,4);

for i = 1:4
  m   =  mass{i};
  x   =  CoM{i}(1); y   =  CoM{i}(2); z   =  CoM{i}(3);
  Ixx =  MoI{i}(1) + m*(y^2+z^2); Iyy =  MoI{i}(2) + m*(x^2+z^2); Izz =  MoI{i}(3) + m*(x^2+y^2);
  Iyz =  PoI{i}(1) - m*y*z; Ixz =  PoI{i}(2) - m*x*z; Ixy =  PoI{i}(3) - m*x*y;

  I(:,:,i)=[(-Ixx+Iyy+Izz)/2, Ixy, Ixz, m*x
             Ixy, (Ixx-Iyy+Izz)/2,Iyz,m*y
             Ixz, Iyz,(Ixx+Iyy-Izz)/2, m*z 
             m*x  m*y m*z m];
end

%% M Matrix

M = sym(zeros(4));
for i = 1:Links
    for j = 1:Links
        for k = max(j,i):Links 
            M(j,i) = M(j,i) + trace(diff(T(:,:,k),q(j))*I(:,:,k)*diff(T(:,:,k)',q(i)));
        end
    end
end
% M=matlabFunction(M);
M=simplify(M);
%% V Matrix
dM_dtheta = sym(zeros(Links^2,Links));
for i = 1:Links
    dM_dtheta((i-1)*Links+1:i*Links,:) = diff(M,q(i));
end

V = simplify(( kron(q_d',eye(Links))-kron(eye(Links),q_d')/2 )*dM_dtheta);
% V= matlabFunction(V);

%% G Matrix
P = sym(0);
G = sym(zeros(Links,1));
for i = 1:Links
   P = P - [g,1]*T(:,:,i)*I(:,:,i)*[0,0,0,1]';
end

for i = 1:Links
      G(i) = diff(P,q(i));
end
% G=matlabFunction(G);
G=simplify(G);

syms t g
syms d(t) th2(t) th3(t) th4(t)
q={d,th2,th3,th4};
I={zeros(3) zeros(3) zeros(3) zeros(3) zeros(3)};
m={0 0 0 0 0};
T={0 0 0 0 0};
td={0 0 0 0 0};
tdd={0 0 0 0 0};
dd={0 0 0 0 0};
ddd={0 0 0 0 0};
P={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
Pc={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
R={zeros(3) zeros(3) zeros(3) zeros(3) zeros(3)};
w={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
wd={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
vd={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
vcd={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
F={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
N={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
f={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
n={zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)};
%mass of links
m0=40.117;
m{1}=22.064;
m{2}=9.805;
m{3}=13.302;
m{4}=10.723;
% center of mass in it's cordination
Pc0=[0;-0.033;0];
Pc{1}=[-0.031;0.165;0.186];
Pc{2}=[-0.005;-0.001;0.082];
Pc{3}=[0.118;0.011;-0.001];
Pc{4}=[0.143;0.015;-0.001];
% moment inertia of links in their cordinations
I0=[13.366 0 0;0 13.428 0;0 0 0.066];
I{1}=[0.134 -0.012 0.001;-0.012 0.165 0;0.001 0 0.186];
I{2}=[0.082 0 0;0 0.053 -0.001;0 -0.001 0.058];
I{3}=[0.048 0.003 0.002;0.003 0.176 0;0.002 0 0.154];
I{4}=[0.018 0.009 -0.001;0.009 0.166 0;-0.001 0 0.164];
% Rotation matrixes
R{1}=[1 0 0;0 1 0;0 0 1];
R{2}=[cos(th2) -sin(th2) 0;0 0 1;-sin(th2) -cos(th2) 0];
R{3}=[-sin(th3) -cos(th3) 0;0 0 -1;cos(th3) -sin(th3) 0];
R{4}=[sin(th4) cos(th4) 0;-cos(th4) sin(th4) 0;0 0 1];
R{5}=[1 0 0;0 1 0;0 0 1];

% each cordination in pervious coordination
P0=[0;0;0];
P{1}=[0;0;d];
P{2}=[0;0.197;0];
P{3}=[0;0;0.1525];
P{4}=[0.290;0;0];
P{5}=[0.3826;0.020;0];
% relative anaular velocities of links
td0=[0;0;0];
td{1}=[0;0;0];
td{2}=[0;0;diff(th2,t)];
td{3}=[0;0;diff(th3,t)];
td{4}=[0;0;diff(th4,t)];


% relative angular accelerations
tdd0=[0;0;0];
tdd{1}=[0;0;0];
tdd{2}=[0;0;diff(th2,t,t)];
tdd{3}=[0;0;diff(th3,t,t)];
tdd{4}=[0;0;diff(th4,t,t)];

% relative Linear velocities of links

% relative Linear velocities of links
dd{1}=[0;0;diff(d,t)];
ddd{1}=[0;0;diff(d,t,t)];

w0=[0;0;0];% initial angular velocity
wd0=[0;0;0];% initial angular acceleration
vd0=[0;g;0];% initial acceleration
f0=[0;0;0];% external force
n0=[0;0;0];% external torques



    
w{1}=R{1}.'*w0+td{1};
wd{1}=R{1}.'*w0;
vd{1}=R{1}.'*(cross(wd0,P{1})+cross(w0,cross(w0,P{1}))+vd0)+cross(2*w{1},dd{1})+ddd{1};
vcd{1}=cross(wd{1},Pc{1})+cross(w{1},cross(w{1},Pc{1}))+vd{1};
F{1}=m{1}*vcd{1};
N{1}=I{1}*wd{1}+cross(w{1},I{1}*w{1});

w{2}=R{2}.'*w{1}+td{2};
wd{2}=R{2}.'*wd{1}+cross(R{2}.'*w{1},td{2})+tdd{2};
vd{2}=R{2}.'*(cross(wd{1},P{2})+cross(w{1},cross(w{1},P{2}))+vd{1});
vcd{2}=cross(wd{2},Pc{2})+cross(w{2},cross(w{2},Pc{2}))+vd{2};
F{2}=m{2}*vcd{2};
N{2}=I{2}*wd{2}+cross(w{2},I{2}*w{2});

w{3}=R{3}.'*w{2}+td{3};
wd{3}=R{3}.'*wd{2}+cross(R{3}.'*w{2},td{3})+tdd{3};
vd{3}=R{3}.'*(cross(wd{2},P{3})+cross(w{2},cross(w{2},P{3}))+vd{2});
vcd{3}=cross(wd{3},Pc{3})+cross(w{3},cross(w{3},Pc{3}))+vd{3};
F{3}=m{3}*vcd{3};
N{3}=I{3}*wd{3}+cross(w{3},I{3}*w{3});

w{4}=R{4}.'*w{3}+td{4};
wd{4}=R{4}.'*wd{3}+cross(R{4}.'*w{3},td{4})+tdd{4};
vd{4}=R{4}.'*(cross(wd{3},P{4})+cross(w{3},cross(w{3},P{4}))+vd{3});
vcd{4}=cross(wd{4},Pc{4})+cross(w{4},cross(w{4},Pc{4}))+vd{4};

F{4}=m{4}*vcd{4};
N{4}=I{4}*wd{4}+cross(w{4},I{4}*w{4});

f{4}=0+F{4};
n{4}=N{4}+0+cross(Pc{4},F{4})+0;
T{4}=n{4}.'*[0;0;1];

f{3}=R{4}*f{4}+F{3};
n{3}=N{3}+R{4}*n{4}+cross(Pc{3},F{3})+cross(P{4},R{4}*f{4});
T{3}=n{3}.'*[0;0;1];

f{2}=R{3}*f{3}+F{2};
n{2}=N{2}+R{3}*n{3}+cross(Pc{2},F{2})+cross(P{3},R{3}*f{3});
T{2}=n{2}.'*[0;0;1];

f{1}=formula(R{2}*f{2}+F{1});
n{1}=N{1}+R{2}*n{2}+cross(Pc{1},F{1})+cross(P{2},R{2}*f{2});
T{1}=n{1}.'*[0;0;1];

FT={f{1}(3) T{2} T{3} T{4}};

for i=1:4
    for j=1:4
        M(i,j)=diff(FT{i},diff(q{j},t,t));
        M(i,j)=simplify(M(i,j));

    end
    G(i,1)=subs(FT{i},{diff(d,t),diff(th2,t),diff(th3,t),diff(th4,t),diff(d,t,t),diff(th2,t,t),diff(th3,t,t),diff(th4,t,t)},{0,0,0,0,0,0,0,0});
    G(i,1)=simplify(G(i,1));
    V(i,1)=subs(FT{i},{diff(d,t,t),diff(th2,t,t),diff(th3,t,t),diff(th4,t,t)},{0,0,0,0})-G(i,1);
    V(i,1)=simplify(V(i,1));
end

save('DynMats_Iterative');

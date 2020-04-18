% Attitude kinematics and dynamics
% the coupling terms do not include RW momentum

function dy = attitudedynamics_Benchmark(t,y)
global MOI_Unc U_RW MA
global DistTorq

wx=y(1);wy=y(2);wz=y(3);
q1=y(4);q2=y(5);q3=y(6);q4=y(7);
h1=y(8);h2=y(9);h3=y(10);h4=y(11);

w=[wx;wy;wz];
wcross=[0,-wz, wy; wz, 0,-wx; -wy, wx, 0]; %wcross=w_cross product
W=[0, wz, -wy, wx; -wz, 0, wx, wy; wy, -wx, 0, wz; -wx, -wy, -wz, 0]; %dq=1/2*W*q
q=[q1;q2;q3;q4];
%u=[ux;uy;uz];

%Total Angular Momentum
H = MOI_Unc*w;

%Attitude Dynamics
dw=inv(MOI_Unc)*(-wcross*H+MA*U_RW+DistTorq);
dq=1/2*W*q;

dy = zeros(7,1);    % a column vector
dy(1) = dw(1);
dy(2) = dw(2);
dy(3) = dw(3);

dy(4) = dq(1);
dy(5) = dq(2);
dy(6) = dq(3);
dy(7) = dq(4);

dy(8) =  -U_RW(1);
dy(9) =  -U_RW(2);
dy(10) = -U_RW(3);
dy(11) = -U_RW(4);
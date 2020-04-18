function dy = Quat_Desire_Propagator(t,y)
global rate_desire

q1_d = y(1);
q2_d = y(2);
q3_d = y(3);
q4_d = y(4);

q_d_Vec = [q1_d, q2_d, q3_d]';
q_d_Scal = q4_d;

q_d_Vec_x = SkewSymmetric(q_d_Vec);

dq_d_Vec = 0.5*(q_d_Vec_x+q_d_Scal*eye(3))*rate_desire;
dq_d_Scal = -0.5*q_d_Vec'*rate_desire;

dy = zeros(4,1);    % a column vector

dy(1) = dq_d_Vec(1);
dy(2) = dq_d_Vec(2);
dy(3) = dq_d_Vec(3);
dy(4) = dq_d_Scal;

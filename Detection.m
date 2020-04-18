function dy = Detection( t, y, omega )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global MOI_Unc MA U_RW_Command

Delta = 5*eye(3);

omega_hat1 = y(1);
omega_hat2 = y(2);
omega_hat3 = y(3);
omega_hat = [omega_hat1; omega_hat2; omega_hat3];

d_omega_hat =  inv(MOI_Unc)* ( - SkewSymmetric(omega_hat)*MOI_Unc*omega_hat ...
                               +  MA*U_RW_Command + Delta*(omega - omega_hat) );

dy = zeros(3,1);                         
dy(1) = d_omega_hat(1);
dy(2) = d_omega_hat(2);
dy(3) = d_omega_hat(3);

end


function dy = Estimation( t, y, additional_para)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global MOI_Unc MA U_RW_Command Est_L Est_G

f_hat = additional_para(1:3)'; 
omega = additional_para(4:6)'; 

omega_hat1 = y(1);
omega_hat2 = y(2);
omega_hat3 = y(3);
omega_hat = [omega_hat1; omega_hat2; omega_hat3];

aux_hat1 = y(4);
aux_hat2 = y(5);
aux_hat3 = y(6);
aux_hat = [aux_hat1; aux_hat2; aux_hat3]; 

d_omega_hat =  inv(MOI_Unc)* ( - SkewSymmetric(omega_hat)*MOI_Unc*omega_hat ...
                               +  MA*U_RW_Command + f_hat + Est_L*(omega-omega_hat) ) ;
d_aux_hat = - Est_G*aux_hat - Est_G*(- SkewSymmetric(omega_hat)*MOI_Unc*omega_hat ... 
                         + MA*U_RW_Command + Est_G*MOI_Unc* omega_hat ) ;                

dy = zeros(6,1);                         
dy(1) = d_omega_hat(1);
dy(2) = d_omega_hat(2);
dy(3) = d_omega_hat(3);
dy(4) = d_aux_hat(1);
dy(5) = d_aux_hat(2);
dy(6) = d_aux_hat(3);

end


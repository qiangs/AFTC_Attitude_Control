function Qe = QuatMultiplication(Qd,Q)

Qd_Vec = Qd(1:3);
Qd_Scal = Qd(4);

Q_Vec = Q(1:3);
Q_Scal = Q(4);

Qd_Vec_x = SkewSymmetric(Qd_Vec);

Qe_Vec = Q_Scal*Qd_Vec+Qd_Scal*Q_Vec+Qd_Vec_x*Q_Vec;
Qe_Scal = Qd_Scal*Q_Scal-Qd_Vec'*Q_Vec;

Qe = [Qe_Vec ; Qe_Scal];
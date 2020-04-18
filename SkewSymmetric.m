function omega_x = SkewSymmetric(omega)

omega1 = omega(1);
omega2 = omega(2);
omega3 = omega(3);

omega_x = [ 0           -omega3         omega2;
            omega3      0               -omega1;
            -omega2     omega1          0];
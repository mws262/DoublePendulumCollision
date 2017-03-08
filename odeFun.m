function zdot = odeFun( t,z,p )

th1dot = z(3);
th2dot = z(4);

% Pick the dynamics based on the contact configuration
switch p.contact
    case 1 % Swings freely
        [th1dotdot,th2dotdot] = accels(p.I1,p.I2,p.d1,p.d2,p.g,p.l1,p.m1,p.m2,z(1),z(2),z(3),z(4));
    case 2 % Link 1 touches
        [th1dotdot,th2dotdot] = accelC1(p.I2,p.d2,p.g,p.m2,z(1),z(2));
    case 3 % Link 2 touches
        [th1dotdot,th2dotdot] = accelC2(p.I1,p.I2,p.d1,p.d2,p.g,p.l1,p.m1,p.m2,p.mu,z(1),z(2),z(3),z(4),p.xcol,p.ycol);
end

zdot = [th1dot;th2dot;th1dotdot;th2dotdot];
end


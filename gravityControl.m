function F = gravityControl(x)
    g = 9.81;
    
    d1 = x(1);
    theta_2 = x(2);
    theta_3 = x(3);

    m2 = 2;
    m3 = 0.6;
    lc2 = 0.5;
    lc3 = 0.5;
    a2 = 1;
    c2 = cos(theta_2);
    c23 = cos(theta_2+theta_3);

    tao2 = m2*g*lc2*c2+m3*g*(a2*c2+lc3*c23);
    tao3 = m3*g*lc3*c23;

    F = [tao2 tao3];
end
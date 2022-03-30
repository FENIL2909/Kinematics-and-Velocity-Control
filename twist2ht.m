function T = twist2ht(S,theta)
format longg
    omega = [S(1) S(2) S(3)];
    v = [S(4) S(5) S(6)]';
    omega_bracket = [0,-omega(3),omega(2); omega(3),0,-omega(1); -omega(2), omega(1),0];
    omega_norm = norm(omega);
    vel_norm = norm(v);
    c = [0 0 0];
    d = [1];
    I = eye(3);
    
    
    if omega_norm >= 0.999999999999
        a = axisangle2rot(omega,theta);
        b = ((I*theta) + ((1-cos(theta))*omega_bracket) + ((theta - sin(theta))*(omega_bracket*omega_bracket)))*v;
        T = [a,b;c,d];
        
    elseif (omega_norm <= 0.0000000000001) && (vel_norm >= 0.999999999999) 
        a = I;
        b = v*theta;
        T = [a,b;c,d];
    end
    
    % If needed, you can calculate a rotation matrix with:
    % R = axisangle2rot(omega,theta);
end
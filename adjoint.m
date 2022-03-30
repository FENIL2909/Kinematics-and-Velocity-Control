function twist_inB = adjoint(twist_inA,T_AB)
    R = [T_AB(1,1),T_AB(1,2),T_AB(1,3); T_AB(2,1), T_AB(2,2), T_AB(2,3); T_AB(3,1), T_AB(3,2), T_AB(3,3)];
    p = [T_AB(1,4), T_AB(2,4), T_AB(3,4)]';
    
    p_bracket = [0,-p(3), p(2); p(3), 0, -p(1); -p(2), p(1), 0];
    K = p_bracket*R;
    O = [0 0 0; 0 0 0; 0 0 0];
    AD_T = [R,O;K,R];
    
    twist_inB = AD_T*twist_inA;
end


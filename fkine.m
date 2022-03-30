function T = fkine(S,M,q)
    SIZE  = size(S);
    T = eye(4);
    
    for i = 1:SIZE(2)
        a = twist2ht(S(:,i),q(i));
        T = T*a;
    end
    T = T*M;
    % If needed, you can convert twists to homogeneous transformation matrices with:
    % twist2ht(S(i),q(i));
end
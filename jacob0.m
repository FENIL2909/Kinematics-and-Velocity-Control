function J = jacob0(S,q) 
SIZE = size(S);
    TRANSFORMATIONS = containers.Map(1,eye(4));
    
    for i = 2:SIZE(2)
        TRANSFORMATIONS(i) = TRANSFORMATIONS(i-1)*twist2ht(S(:,i-1),q(i-1));
    end
    
    
   for j = 1:SIZE(2)
       J(:,j) = adjoint(S(:,j),TRANSFORMATIONS(j));
   end
end
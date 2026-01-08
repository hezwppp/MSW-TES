
% Fit the relationship between MMD and the minimum emissivity
function ychat = Function_MMD(beta0,X)
    
    a0=beta0(1);
    a1=beta0(2);
    a2=beta0(3);
    
    ychat=a0-a1*(X).^a2;

end
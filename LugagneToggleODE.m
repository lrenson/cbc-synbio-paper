function [dx,y] = LugagneToggleODE(t, x, u, klm0, klm, thetaAtc, etaAtc, thetaTet, etaTet, glm, ktm0, ktm, thetaIptg, etaIptg, thetaLac, etaLac, gtm, klp, glp, ktp, gtp, fileArgument)

    %system evolution
    dx=zeros(4,1);

    dx(1)=klm0+klm*HillFunc(x(4)*HillFunc(u(1),thetaAtc,etaAtc),thetaTet,etaTet)-glm*x(1);

    dx(2)=ktm0+ktm*HillFunc(x(3)*HillFunc(u(2),thetaIptg,etaIptg),thetaLac,etaLac)-gtm*x(2);

    dx(3)=klp*x(1)-glp*x(3);

    dx(4)=ktp*x(2)-gtp*x(4);
    
    %outputs TetR LacI
    
    y(1) = x(3);
    
    y(2) = x(4);
end

    function [hill] = HillFunc(x,th,eta)

    hill=1/(1+(x/th)^eta);

end


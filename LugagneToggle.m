    function dx = LugagneToggle(x,u)

    %Parameter Evaluation
    klm0=3.20e-2;

    klm=8.30;

    thetaAtc=11.65;

    etaAtc=2.00;

    thetaTet=30.00;

    etaTet=2.00;

    glm=1.386e-1;

    ktm0=1.19e-1;

    ktm=2.06;

    thetaIptg=9.06e-2;

    etaIptg=2.00;

    thetaLac=31.94;

    etaLac=2.00;

    gtm=1.386e-1;

    klp=9.726e-1;

    glp=1.65e-2;

    ktp=1.170;

    gtp=1.65e-2;


    %system evolution (from the article)
    dx=zeros(4,1);

    dx(1)=klm0+klm*HillFunc(x(4)*HillFunc(u(1),thetaAtc,etaAtc),thetaTet,etaTet)-glm*x(1);

    dx(2)=ktm0+ktm*HillFunc(x(3)*HillFunc(u(2),thetaIptg,etaIptg),thetaLac,etaLac)-gtm*x(2);

    dx(3)=klp*x(1)-glp*x(3);

    dx(4)=ktp*x(2)-gtp*x(4);

    end


    function [hill] = HillFunc(x,th,eta)

    hill=1/(1+(x/th)^eta);

end



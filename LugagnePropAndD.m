function [At,dIPTG,dAtc] = LugagnePropAndD(t,X,Par,Atc,IPTG,extAtc,extIPTG)

    
    At=LugagnePropensities(X,Par,IPTG,Atc);
    
    if extIPTG>IPTG
        
        dIPTG=Par(19)*(extIPTG-IPTG); %IPTG evolution (diffusive term as in the paper)
    else
        dIPTG=Par(21)*(extIPTG-IPTG); %IPTG evolution (diffusive term as in the paper)
    end
    
    if extAtc>Atc
        
        dAtc=Par(20)*(extAtc-Atc);  %Atc Evolution
    else
        dAtc=Par(22)*(extAtc-Atc);  %Atc Evolution
    end
    
    
end

%%
function A = LugagnePropensities(x,Par,IPTG,Atc)
    %Taking the parameters 
    klm0=Par(1);
    klm=Par(2);
    thetaAtc=Par(3);
    etaAtc=Par(4);
    thetaTet=Par(5);
    etaTet=Par(6);
    glm=Par(7);
    ktm0=Par(8);
    ktm=Par(9);   
    thetaIptg=Par(10);
    etaIptg=Par(11);
    thetaLac=Par(12);
    etaLac=Par(13);
    gtm=Par(14);
    klp=Par(15);
    glp=Par(16);
    ktp=Par(17);
    gtp=Par(18);
    
    %Initializing the propensity matrix
    A=zeros(8,1);
    
   
    %Evaluatng the propensity functions (From the article)
    A(1)=glm*x(1);
    A(2)=gtm*x(2);
    A(3)=glp*x(3);
    A(4)=gtp*x(4);
    A(5)=klp*x(1);
    A(6)=ktp*x(2);
    A(7)=klm0+klm*HillFunc(x(4)*HillFunc(Atc,thetaAtc,etaAtc),thetaTet,etaTet);
    A(8)=ktm0+ktm*HillFunc(x(3)*HillFunc(IPTG,thetaIptg,etaIptg),thetaLac,etaLac);

end

    function [hill] = HillFunc(x,th,eta)

    hill=1/(1+(x/th)^eta);

end



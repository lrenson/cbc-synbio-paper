function [tout, xout, Atco, IPTGo] = SDESolver(sim, x0, tspan, Atc0, IPTG0, ext_IPTG, ext_Atc)
    %Sets the initial conditions for x and drugs (internal and external)
    x=x0;
    Atc=Atc0;
    IPTG=IPTG0;
    extAtc=ext_Atc;
    extIPTG=ext_IPTG;
    
    %Control param and time param 
    dt=sim.solvertime;
    t=tspan(1);
    tfin=tspan(2)-2*dt;
    Par=sim.Par;
    sigma=sim.sigma/100;

    S=[-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1;0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0]';
    
    %set the first sample
    xout=x;
    tout=t;
    
    for j=t:dt:tfin
       
        
       [Ak,dIPTG,dAtc]=LugagnePropAndD(t,x,Par,Atc,IPTG,extAtc,extIPTG);
         
        dw=sqrt(dt)*randn(size(S,2),1)*sigma;
        
        x=x+S*Ak*dt+S*diag(sqrt(Ak))*dw;
        
        Atc=Atc+dAtc*dt;
        IPTG=IPTG+dIPTG*dt;
        
        %We update the current time
        time=j+dt;
        
        %Saturate to zero
        for i=1:length(x)
            if(x(i)<0)
               x(i)=0; 
            end
        end
                
        %we set the output variables
        tout=[tout time];
        xout=[xout x];
 
    end
    
    Atco=Atc;
    IPTGo=IPTG;

 
end


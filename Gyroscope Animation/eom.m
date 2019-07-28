function [eoms] = eom()
    % Init    
    syms t a(t) b(t) c(t) d(t) D L R MGx MGy FGx FGy FGz FOx FOy FOz adot bdot cdot ddot...
        addot bddot cddot dddot
    R1to0=rMat('z','p','a(t)');
    R0to1=R1to0.';
    R2to1=rMat('x','p','b(t)');
    R1to2=R2to1.';
    R3to2=rMat('z','p','c(t)');
    R2to3=R3to2.';
    R4to3=rMat('z','p','d(t)');
    R3to4=R4to3.';
    X=zeros(1,8);
    Xdot=zeros(1,8);
    r=2.5/2; %frame radius
    Ro=70/2; %frame width
    Ri=57/2; %rotor width
    H=94; %OC height
    h=6.5; %rotor thickness
    L=(H/2)+Ro; %OG distance
    mr=45; %measured mass of rotor
    rv=h*pi()*((Ri)^2); %rotor volume
    rho=mr/rv; 
    g=9.81;

    %A1
    w32_3=[0;0;diff(c(t))];
    w21_2=[diff(b(t));0;0];
    w21_3=R2to3*w21_2;
    w10_1=[0;0;diff(a(t))];
    w10_3=R2to3*R1to2*w10_1;
    w3_3=w32_3+w21_3+w10_3;

    w3_0=R1to0*R2to1*R3to2*w3_3;
    w3_3dot=diff(w3_3);
    w43_3=[0;0;diff(d(t))];
    w4_3=w43_3+w3_3;
    w4_3dot=diff(w4_3)+cross(w3_3,w4_3);

    rOG_2=[0;0;L];
    rOG_0=R1to0*R2to1*rOG_2;
    vOG_0=diff(rOG_0);
    rAB_3=[D;0;0];
    rGP_4=[R;0;0];

    rGP_3=R4to3*rGP_4;
    vGP_3=diff(rGP_3) +cross(w3_3,rGP_3);
    aGP_3=diff(vGP_3)+cross(w3_3,vGP_3);
    
    %A2
    Ir_G_3=mr*[(3*(Ri^2)+h^2)/12 0 0;0 (3*(Ri^2)+h^2)/12 0;0 0 (Ri^2)/2];
    
    mar=H*rho*pi()*r^2;
    mbr=2*rho*Ro*(pi()^2)*r^2;
    mcr=mbr;
    mf=23;
    mt=mar+mbr+mcr;
    ma=mf*(mar/mt);
    mb=mf*(mbr/mt);
    mc=mf*(mcr/mt);
    
    Ia_G_3=ma*[(3*(r^2)+H^2)/12 0 0;0 (3*(r^2)+H^2)/12 0;0 0 (r^2)/2];
    Ib_G_3=mb*[(5/8)*r^2 + (1/2)*Ro^2 0 0;0 (3/4)*r^2 + Ro^2 0;0 0 ...
        (5/8)*r^2 + (1/2)*Ro^2];
    Ic_G_3=mc*[(5/8)*r^2 + (1/2)*Ro^2 0 0;0 (3/4)*r^2 + Ro^2 0;0 0 ...
        (5/8)*r^2 + (1/2)*Ro^2];
    
    If_G_3=Ia_G_3 + Ib_G_3 + Ic_G_3;
    
    rOG_3=[0;0;L];
    delta=deltaMat(rOG_3);
    If_O_3=If_G_3 + mf*delta;

    rOG_3ddot(1)=(-L*diff(b(t),t,t)*sin(c(t))) + (L*diff(a(t),t,t)*sin(b(t))...
        *cos(c(t))) + (2*L*diff(a(t),t)*diff(b(t),t)*cos(b(t))*cos(c(t)))...
        +(L*cos(b(t))*sin(b(t))*sin(c(t))*(diff(a(t),t)^2));
    rOG_3ddot(2)=...
        (-L*diff(b(t),t,t)*cos(c(t)))-(2*L*diff(a(t),t)*diff(b(t),t)...
        *cos(b(t))*sin(c(t))) + (L*cos(b(t))*sin(b(t))*cos(c(t))*(diff(a(t),t)^2));
    rOG_3ddot(3)=...
        -(L*(diff(b(t),t)^2))-(L*(diff(a(t),t)^2)*sin(b(t))*sin(b(t)));
    
    %Rotor EOM
    Frotor=[FGx-(490*sin(b(t))*sin(c(t)));FGy-(490*sin(b(t))*cos(c(t)));...
        FGz-(490*cos(b(t)))];
    FGx=solve(Frotor(1)==(mr*rOG_3ddot(1)),FGx);
    FGy=solve(Frotor(2)==(mr*rOG_3ddot(2)),FGy);
    FGz=solve(Frotor(3)==(mr*rOG_3ddot(3)),FGz);
    
    Mrotor=[MGx;MGy;0];
    %w4=wrotor
    hrotor_G_3=Ir_G_3*w4_3;
    hrotor_G_3dot=diff(hrotor_G_3)+cross(w3_3,hrotor_G_3);
    MGx=solve(Mrotor(1)==hrotor_G_3dot(1),MGx);
    MGy=solve(Mrotor(2)==hrotor_G_3dot(2),MGy);
    
    %Frame EOM
    Ff=[FOx-FGx-(mf*g*sin(b(t))*sin(c(t)));...
        FOy-FGy-(mf*g*sin(b(t))*cos(c(t)));...
        FOz-FGz-(mf*g*cos(b(t)))];
    FOx=solve(Ff(1)==(mf*rOG_3ddot(1)),FOx);
    FOy=solve(Ff(2)==(mf*rOG_3ddot(2)),FOy);
    FOz=solve(Ff(3)==(mf*rOG_3ddot(3)),FOz);
    
    hf_O_3=If_O_3*w3_3;
    hf_O_3dot=diff(hf_O_3) + cross(w3_3,hf_O_3);
    Gf_1=[0;0;-mf*g];
    Gf_3=R2to3*R1to2*Gf_1;
    Fg_3=[FGx;FGy;FGz];
    Mrotor=[MGx;MGy;0];
    Mf_O_3=-Mrotor+cross(rOG_3,Gf_3 -Fg_3);
    
    %Set RHS = 0 for double derivate eqns
    EOM=[hrotor_G_3dot(3), Mf_O_3(1)-hf_O_3dot(1), Mf_O_3(2)-hf_O_3dot(2), Mf_O_3(3)-hf_O_3dot(3)];
    %Substituting out diff() terms as they can't be solved for
    EOM=subs(EOM, [diff(a(t),t) diff(b(t),t) diff(c(t),t) diff(d(t),t)...
        diff(a(t),t,t) diff(b(t),t,t) diff(c(t),t,t) diff(d(t),t,t)], [adot bdot cdot ddot...
        addot bddot cddot dddot]);
    %dderivs is a matrix containing final EOM
    dderivs=solve([EOM(2)==0 EOM(3)==0 EOM(4)==0 EOM(1)==0],[addot...
        bddot cddot dddot]);
    eoms=[dderivs.addot ; dderivs.bddot ; dderivs.cddot ; dderivs.dddot];
    
    
end
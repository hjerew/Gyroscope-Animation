clear all
syms t x1 x2 x3 x4 x5 x6 x7 x8 a(t) b(t) c(t) d(t) adot bdot cdot ddot

R1to0=rMat('z','p','a(t)');
R0to1=R1to0.';
R2to1=rMat('x','p','b(t)');
R1to2=R2to1.';
R3to2=rMat('z','p','c(t)');
R2to3=R3to2.';  
R4to3=rMat('z','p','d(t)');
R3to4=R4to3.';

eoms=eom();
eoms=subs(eoms, [a(t) b(t) c(t) d(t) adot bdot cdot ddot], [x1 x2 x3 x4 x5 x6 x7 x8]);


x_init = [0.01;0.01;0.01;0;20;0.01;0.01;20];  %[a;b;c;d;a_dot;b_dot;c_dot;d_dot]  
t_span = [0 2000];                                 % start and end times
options = odeset('RelTol',1e-8,'AbsTol',1e-8'); % solver options
X_sol = ode45(@stateSpace,t_span,x_init, options)         % solve the eoms
VIDEO = 1;

% ------------------------------------------------------------------------
%Evaluation the solution

dt = 0.03;                                  % set time step                        
t = t_span(1):dt:t_span(2);                % creat vector of time
A = deval(X_sol,t);                           
display(length(t))

% ------------------------------------------------------------------------
% Plotting the states
plot(t,A)
xlabel('Time')
ylabel('States')
h = legend('$\alpha$','$\beta$','$\gamma$','$\delta$','$\dot{\alpha}$',...
    '$\dot{\beta}$','$\dot{\gamma}$','$\dot{\delta}$');
set(h,'Interpreter','latex')

% ------------------------------------------------------------------------
% Gyroscope Creation
 
%  FRAME (Torque)
t = 0:2*pi/100:2*pi;
x1 = 16*sin(t);
y1 = 16*cos(t);
z1 = linspace(15,15,length(x1));
 
%  FRAME (Torque)
x2 = 16*sin(t);
y2 =linspace(0,0,length(x1));
z2 = 18*cos(t)+15;
 
%  FRAME a (Cylinder)
x3 = linspace(0,0,length(x1));
y3 = linspace(0,0,length(x1));
z3 = linspace(-10,40,length(x1));

%[x3,y3,z3]=cylinder([0 0 0 0 14 14 0 0 0 0])
%  ROTOR
[x4,y4,z4] = cylinder([0 0 0 0 14 14 0 0 0 0]);
z4 = 10*z4+10;
 
% ------------------------------------------------------------------------
% Setup Video

if VIDEO
    fps =30;% 1/dt;
    Video = VideoWriter('Group 69','MPEG-4');
    Video.FrameRate = fps;
    open(Video);
end
 

% ------------------------------------------------------------------------
%Creating Animation

Gyroscope = figure;
hold on
for i = 1:length(t)
    cla
    
    alpha = A(1,i);  %a=alpha
    beta = A(2,i);  %b=Beta
    gamma= A(3,i);   %c=gamma
    delta = A(4,i);  %d=delta
    
    %rotation matrices
     R01 = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
    %R01=subs(R0to1, a, alpha);
     R12 = [1 0 0; 0 cos(beta) sin(beta); 0 -sin(beta) cos(beta)];
    %R12=subs(R1to2, b, beta);
     R23 = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
    %R23=subs(R2to3, c, gamma);
     R34 = [cos(delta) sin(delta) 0; -sin(delta) cos(delta) 0; 0 0 1];
    %R34=subs(R3to4, d, delta);
    
    %rotated gyroscope
    [x1_rotated,y1_rotated,z1_rotated] = rotate(x1,y1,z1,R01*R12*R23);
    [x2_rotated,y2_rotated,z2_rotated] = rotate(x2,y2,z2,R01*R12*R23);
    [x3_rotated,y3_rotated,z3_rotated] = rotate(x3,y3,z3,R01*R12*R23);
    [x4_rotated,y4_rotated,z4_rotated] = rotate(x4,y4,z4,R01*R12*R23*R34);
    
    plot3(x1_rotated,y1_rotated,z1_rotated,'linewidth',3.5)
    plot3(x2_rotated,y2_rotated,z2_rotated,'linewidth',3.5)
    plot3(x3_rotated,y3_rotated,z3_rotated,'linewidth',3.5)
   surf(x4_rotated,y4_rotated,z4_rotated,x4,'linewidth', 2)
    
    axis square
    view(45,25)
    axis(1*[-20 20 -20 20 -10 40])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    if VIDEO
        writeVideo(Video,getframe(Gyroscope));
    else
        pause(dt)
    end
end
 
if VIDEO
close(Video)
end
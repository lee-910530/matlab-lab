clc;
close all;
clear all;
%..........
%.......data........
h = 0.05;
fin = 3000;
tn = 0:h:fin

G_pat1(1) = 70;
X_pat1(1) = 0;
I_pat1(1) = 7;
G_pat2(1) = 70;
X_pat2(1) = 0;
I_pat2(1) = 7;
G_pat3(1) = 70;
X_pat3(1) = 0;
I_pat3(1) = 7;

G_pat1_d = zeros(1,fin/h+1);
U1 = zeros(1,fin/h+1);
G_pat2_d = zeros(1,fin/h+1);
U2 = zeros(1,fin/h+1);
G_pat3_d = zeros(1,fin/h+1);
U3 = zeros(1,fin/h+1);


x = fin/h;
%..........glucose concentration-input1...........
x1 = 0:0.02:400;
NB = trapezoid_up_wave(x1,50,70);
NS = triangle_wave(x1,50,70,75);
NOM = triangle_wave(x1,70,75,80);
PS = triangle_wave(x1,75,80,150);
PM = triangle_wave(x1,75,150,225);
PB = triangle_wave(x1,150,225,300);
PL = trapezoid_down_wave(x1,225,300);

% figure('name','glucose concentration-input1');
% plot(x1,NB,"b",x1,NS,"g",x1,NOM,"r",x1,PS,"c",x1,PL,"k");
% hold on;
% plot(x1,PM,"Color","#E020F0")
% hold on;
% plot(x1,PB,"Color","#E09010")
% xlim([40,400])
% title('glucose concentration-input1')
% xlabel("glucose concentration(mg/dl)")

%..........glucose deviation-input2...........
x2 = -20:0.002:20;
neg = triangle_wave(x2,-20,-20,0);
zero = triangle_wave(x2,-1,0,1);
pos = triangle_wave(x2,0,20,20);

% figure('name','glucose deviation-input2');
% plot(x2,neg,"b",x2,zero,"g",x2,pos,"r");
% title('glucose deviation-input2')
% xlabel("glucose deviation(mg/dl)")

%..........insulin infusion rate-output...........
y = -1:0.0005:9;
NB_o = triangle_wave(y,-1,-0.5,-0.2);
NS_o = triangle_wave(y,-0.4,-0.2,0.2);
Z_o = triangle_wave(y,-0.2,0.2,0.5);
PS_o = triangle_wave(y,0.2,0.5,2);
PM_o = triangle_wave(y,0.6,2,4);
PB_o = triangle_wave(y,2,4,6);
PL_o = triangle_wave(y,4,6,8);

% figure('name','insulin infusion rate-output');
% plot(y,NB_o,"b",y,NS_o,"g",y,Z_o,"r",y,PS_o,"c",y,PL_o,"k");
% hold on;
% plot(y,PM_o,"Color","#E020F0")
% hold on;
% plot(y,PB_o,"Color","#E09010")
% title("insulin infusion rate-output")
% xlabel("insulin infusion rate(microU/min/mg)")

%.............patient3................
num = 1;
for i = tn
    G_pat1_d(num) = G_p1(i,G_pat1(num),X_pat1(num),I_pat1(num));
    G_pat2_d(num) = G_p2(i,G_pat2(num),X_pat2(num),I_pat2(num));
    G_pat3_d(num) = G_p3(i,G_pat3(num),X_pat3(num),I_pat3(num));

    q1 = G_pat1(num);
    q2 = G_pat1_d(num);

    [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10] = MU(q1,q2);
        
    w1 = min(mu7,mu8);  %PL
    w2 = min(mu7,mu9); %PL
    w3 = min(mu7,mu10);  %PL
    w4 = min(mu6,mu8);  %PB
    w5 = min(mu6,mu9); %PB
    w6 = min(mu6,mu10);  %PL
    w7 = min(mu5,mu8);  %PM
    w8 = min(mu5,mu9); %PM
    w9 = min(mu5,mu10);  %PM
    w10 = min(mu4,mu8); %PS
    w11 = min(mu4,mu9);%PS
    w12 = min(mu4,mu10); %PS
    w13 = min(mu3,mu8);%Z
    w14 = min(mu3,mu9);%Z
    w15 = min(mu3,mu10); %Z
    w16 = min(mu2,mu8);  %NB
    w17 = min(mu2,mu9); %NS
    w18 = min(mu2,mu10);  %NS
    w19 = min(mu1,mu8);  %NB
    w20 = min(mu1,mu9); %NB
    w21 = min(mu1,mu10);  %NB

%     NB_o = -0.5;
%     NS_o = -0.2;
%     Z_o = 0.1;
%     PS_o = 0.5;
%     PM_o = 2;
%     PB_o = 4;
%     PL_o = 6;
    NB_o = -0.5;
    NS_o = -0.2;
    Z_o = 0.1;
    PS_o = 0.5;
    PM_o = 2;
    PB_o = 4;
    PL_o = 6;
    
    nn = (w1*PL_o+w2*PL_o+w3*PL_o+w4*PB_o+w5*PB_o+w6*PL_o+w7*PM_o+w8*PM_o+w9*PM_o+w10*PS_o ...
        +w11*PS_o+w12*PS_o+w13*Z_o+w14*Z_o+w15*Z_o+w16*NB_o+w17*NS_o+w18*NS_o+w19*NB_o+w20*NB_o+w21*NB_o)...
        /(w1+w2+w3+w4+w5+w6+w7+w8+w9+w10+w11+w12+w13+w14+w15+w16+w17+w18+w19+w20+w21);

    U1(num) = sat(nn);

    k1 = G_p1(i,G_pat1(num),X_pat1(num),I_pat1(num));
    l1 = X_p1(i,G_pat1(num),X_pat1(num),I_pat1(num));
    j1 = I_p1(i,G_pat1(num),X_pat1(num),I_pat1(num),U1(num));

    k2 = G_p1(i + h/2,G_pat1(num) + h*k1/2,X_pat1(num) + h*l1/2,I_pat1(num) + h*j1/2);
    l2 = X_p1(i + h/2,G_pat1(num) + h*k1/2,X_pat1(num) + h*l1/2,I_pat1(num) + h*j1/2);
    j2 = I_p1(i + h/2,G_pat1(num) + h*k1/2,X_pat1(num) + h*l1/2,I_pat1(num) + h*j1/2,U1(num));

    k3 = G_p1(i + h/2,G_pat1(num) + h*k2/2,X_pat1(num) + h*l2/2,I_pat1(num) + h*j2/2);
    l3 = X_p1(i + h/2,G_pat1(num) + h*k2/2,X_pat1(num) + h*l2/2,I_pat1(num) + h*j2/2);
    j3 = I_p1(i + h/2,G_pat1(num) + h*k2/2,X_pat1(num) + h*l2/2,I_pat1(num) + h*j2/2,U1(num));
    
    k4 = G_p1(i + h,G_pat1(num) + h*k3,X_pat1(num) + h*l3,I_pat1(num) + h*j3);
    l4 = X_p1(i + h,G_pat1(num) + h*k3,X_pat1(num) + h*l3,I_pat1(num) + h*j3);
    j4 = I_p1(i + h,G_pat1(num) + h*k3,X_pat1(num) + h*l3,I_pat1(num) + h*j3,U1(num));
    
    if num < x +1
        G_pat1(num+1) = G_pat1(num) + h*(k1+2*k2+2*k3+k4)/6;
        X_pat1(num+1) = X_pat1(num) + h*(l1+2*l2+2*l3+l4)/6;
        I_pat1(num+1) = I_pat1(num) + h*(j1+2*j2+2*j3+j4)/6;
    end

    q1 = G_pat2(num);
    q2 = G_pat2_d(num);

    [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10] = MU(q1,q2);
        
    w1 = min(mu7,mu8);  %PL
    w2 = min(mu7,mu9); %PL
    w3 = min(mu7,mu10);  %PL
    w4 = min(mu6,mu8);  %PB
    w5 = min(mu6,mu9); %PB
    w6 = min(mu6,mu10);  %PL
    w7 = min(mu5,mu8);  %PM
    w8 = min(mu5,mu9); %PM
    w9 = min(mu5,mu10);  %PM
    w10 = min(mu4,mu8); %PS
    w11 = min(mu4,mu9);%PS
    w12 = min(mu4,mu10); %PS
    w13 = min(mu3,mu8);%Z
    w14 = min(mu3,mu9);%Z
    w15 = min(mu3,mu10); %Z
    w16 = min(mu2,mu8);  %NB
    w17 = min(mu2,mu9); %NS
    w18 = min(mu2,mu10);  %NS
    w19 = min(mu1,mu8);  %NB
    w20 = min(mu1,mu9); %NB
    w21 = min(mu1,mu10);  %NB

    
    nn = (w1*PL_o+w2*PL_o+w3*PL_o+w4*PB_o+w5*PB_o+w6*PL_o+w7*PM_o+w8*PM_o+w9*PM_o+w10*PS_o ...
        +w11*PS_o+w12*PS_o+w13*Z_o+w14*Z_o+w15*Z_o+w16*NB_o+w17*NS_o+w18*NS_o+w19*NB_o+w20*NB_o+w21*NB_o)...
        /(w1+w2+w3+w4+w5+w6+w7+w8+w9+w10+w11+w12+w13+w14+w15+w16+w17+w18+w19+w20+w21);

    U2(num) = sat(nn);

    k1 = G_p2(i,G_pat2(num),X_pat2(num),I_pat2(num));
    l1 = X_p2(i,G_pat2(num),X_pat2(num),I_pat2(num));
    j1 = I_p2(i,G_pat2(num),X_pat2(num),I_pat2(num),U2(num));

    k2 = G_p2(i + h/2,G_pat2(num) + h*k1/2,X_pat2(num) + h*l1/2,I_pat2(num) + h*j1/2);
    l2 = X_p2(i + h/2,G_pat2(num) + h*k1/2,X_pat2(num) + h*l1/2,I_pat2(num) + h*j1/2);
    j2 = I_p2(i + h/2,G_pat2(num) + h*k1/2,X_pat2(num) + h*l1/2,I_pat2(num) + h*j1/2,U2(num));

    k3 = G_p2(i + h/2,G_pat2(num) + h*k2/2,X_pat2(num) + h*l2/2,I_pat2(num) + h*j2/2);
    l3 = X_p2(i + h/2,G_pat2(num) + h*k2/2,X_pat2(num) + h*l2/2,I_pat2(num) + h*j2/2);
    j3 = I_p2(i + h/2,G_pat2(num) + h*k2/2,X_pat2(num) + h*l2/2,I_pat2(num) + h*j2/2,U2(num));
    
    k4 = G_p2(i + h,G_pat2(num) + h*k3,X_pat2(num) + h*l3,I_pat2(num) + h*j3);
    l4 = X_p2(i + h,G_pat2(num) + h*k3,X_pat2(num) + h*l3,I_pat2(num) + h*j3);
    j4 = I_p2(i + h,G_pat2(num) + h*k3,X_pat2(num) + h*l3,I_pat2(num) + h*j3,U2(num));
    
    if num < x +1
        G_pat2(num+1) = G_pat2(num) + h*(k1+2*k2+2*k3+k4)/6;
        X_pat2(num+1) = X_pat2(num) + h*(l1+2*l2+2*l3+l4)/6;
        I_pat2(num+1) = I_pat2(num) + h*(j1+2*j2+2*j3+j4)/6;
    end

    q1 = G_pat3(num);
    q2 = G_pat3_d(num);

    [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10] = MU(q1,q2);
        
    w1 = min(mu7,mu8);  %PL
    w2 = min(mu7,mu9); %PL
    w3 = min(mu7,mu10);  %PL
    w4 = min(mu6,mu8);  %PB
    w5 = min(mu6,mu9); %PB
    w6 = min(mu6,mu10);  %PL
    w7 = min(mu5,mu8);  %PM
    w8 = min(mu5,mu9); %PM
    w9 = min(mu5,mu10);  %PM
    w10 = min(mu4,mu8); %PS
    w11 = min(mu4,mu9);%PS
    w12 = min(mu4,mu10); %PS
    w13 = min(mu3,mu8);%Z
    w14 = min(mu3,mu9);%Z
    w15 = min(mu3,mu10); %Z
    w16 = min(mu2,mu8);  %NB
    w17 = min(mu2,mu9); %NS
    w18 = min(mu2,mu10);  %NS
    w19 = min(mu1,mu8);  %NB
    w20 = min(mu1,mu9); %NB
    w21 = min(mu1,mu10);  %NB

    
    nn = (w1*PL_o+w2*PL_o+w3*PL_o+w4*PB_o+w5*PB_o+w6*PL_o+w7*PM_o+w8*PM_o+w9*PM_o+w10*PS_o ...
        +w11*PS_o+w12*PS_o+w13*Z_o+w14*Z_o+w15*Z_o+w16*NB_o+w17*NS_o+w18*NS_o+w19*NB_o+w20*NB_o+w21*NB_o)...
        /(w1+w2+w3+w4+w5+w6+w7+w8+w9+w10+w11+w12+w13+w14+w15+w16+w17+w18+w19+w20+w21);

    U3(num) = sat(nn);

    k1 = G_p3(i,G_pat3(num),X_pat3(num),I_pat3(num));
    l1 = X_p3(i,G_pat3(num),X_pat3(num),I_pat3(num));
    j1 = I_p3(i,G_pat3(num),X_pat3(num),I_pat3(num),U3(num));

    k2 = G_p3(i + h/2,G_pat3(num) + h*k1/2,X_pat3(num) + h*l1/2,I_pat3(num) + h*j1/2);
    l2 = X_p3(i + h/2,G_pat3(num) + h*k1/2,X_pat3(num) + h*l1/2,I_pat3(num) + h*j1/2);
    j2 = I_p3(i + h/2,G_pat3(num) + h*k1/2,X_pat3(num) + h*l1/2,I_pat3(num) + h*j1/2,U3(num));

    k3 = G_p3(i + h/2,G_pat3(num) + h*k2/2,X_pat3(num) + h*l2/2,I_pat3(num) + h*j2/2);
    l3 = X_p3(i + h/2,G_pat3(num) + h*k2/2,X_pat3(num) + h*l2/2,I_pat3(num) + h*j2/2);
    j3 = I_p3(i + h/2,G_pat3(num) + h*k2/2,X_pat3(num) + h*l2/2,I_pat3(num) + h*j2/2,U3(num));
    
    k4 = G_p3(i + h,G_pat3(num) + h*k3,X_pat3(num) + h*l3,I_pat3(num) + h*j3);
    l4 = X_p3(i + h,G_pat3(num) + h*k3,X_pat3(num) + h*l3,I_pat3(num) + h*j3);
    j4 = I_p3(i + h,G_pat3(num) + h*k3,X_pat3(num) + h*l3,I_pat3(num) + h*j3,U3(num));
    
    if num < x +1
        G_pat3(num+1) = G_pat3(num) + h*(k1+2*k2+2*k3+k4)/6;
        X_pat3(num+1) = X_pat3(num) + h*(l1+2*l2+2*l3+l4)/6;
        I_pat3(num+1) = I_pat3(num) + h*(j1+2*j2+2*j3+j4)/6;
    end
       num = num +1;
end

figure(6)
subplot(3,1,3)
hold on
plot(tn,U1,'b')
plot(tn,U2,'--r')
plot(tn,U3,'--g')
yline(0)
title('u(t)')
ylabel('(MicroU)/mg/min2')
xlabel('time(min)')

subplot(3,1,2)
hold on
plot(tn,I_pat1,'b')
plot(tn,I_pat2,'--r')
plot(tn,I_pat3,'--g')
yline(7)
title('insulin')
xlabel('time(min)')
ylabel('(MicroU/ml)')

subplot(3,1,1)
hold on
plot(tn,G_pat1,'b')
plot(tn,G_pat2,'--r')
plot(tn,G_pat3,'--g')
yline(70)
title('G(t)')
ylabel('(mg/dl)')
xlabel('time(min)')

function g = G_n(t,p1,Gb,G,X,I)
    w = t-500;
    w(w<0) = inf;
    g = -p1*(G-Gb)-X*G + 0.3*exp(-0.02*t)+ 0.3*exp(-0.02*w);
end
  
function x = X_n(t,p2,p3,Ib,G,X,I)
    x = -p2*X+p3*(I-Ib);
end

function i = I_n(t,n,Ib,r,h,G,X,I,u)
    x = max(G-h,0);
    i = -n*(I-Ib)+r*x*t+u;
end


%..........patient1..........
function g_p1 = G_p1(t,G,X,I)
    g_p1 = G_n(t,0,70,G,X,I);
end

function x_p1 = X_p1(t,G,X,I)
    x_p1 = X_n(t,0.02,5.3*10^-6,7,G,X,I);

end

function i_p1 = I_p1(t,G,X,I,u)
 %   i_p1 = I_n(t,0.3,7,0.005,78,G,X,I,u);
     i_p1 = I_n(t,0.3,7,0,78,G,X,I,u);
end

%..........patient2..........
function g_p2 = G_p2(t,G,X,I)
    g_p2 = G_n(t,0,70,G,X,I);
end

function x_p2 = X_p2(t,G,X,I)
    x_p2 = X_n(t,0.0072,2.16*10^-6,7,G,X,I);

end

function i_p2 = I_p2(t,G,X,I,u)
%    i_p2 = I_n(t,0.2456,7,0.0038,77.5783,G,X,I,u);
     i_p2 = I_n(t,0.2456,7,0,77.5783,G,X,I,u);
end

%..........patient3..........
function g_p3 = G_p3(t,G,X,I)
    g_p3 = G_n(t,0,70,G,X,I);
end

function x_p3 = X_p3(t,G,X,I)
    x_p3 = X_n(t,0.0142,9.94*10^-6,7,G,X,I);

end

function i_p3 = I_p3(t,G,X,I,u)
 %   i_p3 = I_n(t,0.2814,7,0.0046,82.9370,G,X,I,u);
     i_p3 = I_n(t,0.2814,7,0,82.9370,G,X,I,u);
end





%..........function...........
function trap = trapezoid_up_wave(t,start,fin)
    num = 0;
    for i = t
        num = num + 1;
        if i < start
            trap(num) = 1;
        elseif i >= fin
            trap(num) = 0;
        elseif i >= start && i < fin
            trap(num) = 1/(start - fin) * (i - fin);
        end
    end
end

function trap = trapezoid_down_wave(t,start,fin)
    num = 0;
    for i = t
        num = num + 1;
        if i < start
            trap(num) = 0;
        elseif i >= fin
            trap(num) = 1;
        elseif i >= start && i < fin
            trap(num) = 1/(fin - start) * (i - start);
        end
    end
end

function tri = triangle_wave(t,start,middle,fin)
    num = 0;
    for i = t
        num = num + 1;
        if i < start
            tri(num) = 0;
        elseif i >= fin
            tri(num) = 0;
        elseif i >= start && i < middle
            tri(num) = 1/(middle-start)*(i-start);
        else 
            tri(num) = 1/(middle-fin)*(i-fin);
        end
    end
end

function mem = tri_member(i,start,middle,fin)
    if i < start
        mem = 0;
    elseif i >= fin
        mem = 0;
    elseif i >= start && i < middle
        mem = 1/(middle-start)*(i-start);
    else 
        mem = 1/(middle-fin)*(i-fin);
    end
end

function mem1 = trape_down_member(i,start,fin)
    if i < start
        mem1 = 0;
    elseif i >= fin
        mem1 = 1;
    elseif i >= start && i < fin
        mem1 = 1/(fin - start) * (i - start);
    end
end

function mem2 = trape_up_member(i,start,fin)
    if i < start
        mem2 = 1;
    elseif i >= fin
        mem2 = 0;
    elseif i >= start && i < fin
        mem2 = 1/(start - fin) * (i - fin);
    end
end

function y=sat(x);
% sat is the saturation function with unit limits and unit slope.
if x>1
    y=1;
elseif x<0 
    y=0;
else 
    y=x;
end
end

function [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,mu9,mu10] = MU(q1,q2)
    mu1 = trape_up_member(q1,40,60);
    mu2 = tri_member(q1,50,60,75);
    mu3 = tri_member(q1,60,75,80);
    mu4 = tri_member(q1,75,78,100);
    mu5 = tri_member(q1,75,150,225);
    mu6 = tri_member(q1,150,225,300);
    mu7 = trape_down_member(q1,225,300);

    mu8 = tri_member(q2,-20,-20,0);
    mu9 = tri_member(q2,-5,0,5);
    mu10 = tri_member(q2,0,20,20);

%     mu1 = trape_up_member(q1,40,60);
%     mu2 = tri_member(q1,50,55,90);
%     mu3 = tri_member(q1,70,75,80);
%     mu4 = tri_member(q1,60,80,100);
%     mu5 = tri_member(q1,75,150,225);
%     mu6 = tri_member(q1,200,250,300);
%     mu7 = trape_down_member(q1,225,400);
% 
%     mu8 = tri_member(q2,-20,-20,0);
%     mu9 = tri_member(q2,-1,0,1);
%     mu10 = tri_member(q2,0,20,20);


end

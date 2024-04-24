clc;
clf;
clear all;
%..........
%.......data........
h = 0.05;
fin = 1000;
tn = 0:h:fin

G_pat3(1) = 220;
X_pat3(1) = 0;
I_pat3(1) = 50;
I3 = zeros(1,fin/h+1);
X3 = zeros(1,fin/h+1);
G3 = zeros(1,fin/h+1);

G_normal(1) = 291.2;
X_normal(1) = 0;
I_normal(1) = 364.8;
In = zeros(1,fin/h+1);
Xn = zeros(1,fin/h+1);
Gn = zeros(1,fin/h+1);

%.......runge kutta patient........
x = fin/h;

num = 1;
for i = tn
    k1 = G_p3(tn(num),G_pat3(num),X_pat3(num),I_pat3(num));
    l1 = X_p3(tn(num),G_pat3(num),X_pat3(num),I_pat3(num));
    j1 = I_p3(tn(num),G_pat3(num),X_pat3(num),I_pat3(num),0);

    k2 = G_p3(tn(num) + h/2,G_pat3(num) + h*k1/2,X_pat3(num) + h*l1/2,I_pat3(num) + h*j1/2);
    l2 = X_p3(tn(num) + h/2,G_pat3(num) + h*k1/2,X_pat3(num) + h*l1/2,I_pat3(num) + h*j1/2);
    j2 = I_p3(tn(num) + h/2,G_pat3(num) + h*k1/2,X_pat3(num) + h*l1/2,I_pat3(num) + h*j1/2,0);

    k3 = G_p3(tn(num) + h/2,G_pat3(num) + h*k2/2,X_pat3(num) + h*l2/2,I_pat3(num) + h*j2/2);
    l3 = X_p3(tn(num) + h/2,G_pat3(num) + h*k2/2,X_pat3(num) + h*l2/2,I_pat3(num) + h*j2/2);
    j3 = I_p3(tn(num) + h/2,G_pat3(num) + h*k2/2,X_pat3(num) + h*l2/2,I_pat3(num) + h*j2/2,0);
    
    k4 = G_p3(tn(num) + h,G_pat3(num) + h*k3,X_pat3(num) + h*l3,I_pat3(num) + h*j3);
    l4 = X_p3(tn(num) + h,G_pat3(num) + h*k3,X_pat3(num) + h*l3,I_pat3(num) + h*j3);
    j4 = I_p3(tn(num) + h,G_pat3(num) + h*k3,X_pat3(num) + h*l3,I_pat3(num) + h*j3,0);

    if num <= x
        G_pat3(num+1) = G_pat3(num) + h*(k1+2*k2+2*k3+k4)/6;
        X_pat3(num+1) = X_pat3(num) + h*(l1+2*l2+2*l3+l4)/6;
        I_pat3(num+1) = I_pat3(num) + h*(j1+2*j2+2*j3+j4)/6;
    
        num = num + 1;
    end
end

%.......runge kutta normal........
num = 1;
for i = tn
    k1 = G_no(tn(num),G_normal(num),X_normal(num),I_normal(num));
    l1 = X_no(tn(num),G_normal(num),X_normal(num),I_normal(num));
    j1 = I_no(tn(num),G_normal(num),X_normal(num),I_normal(num),0);

    k2 = G_no(tn(num) + h/2,G_normal(num) + h*k1/2,X_normal(num) + h*l1/2,I_normal(num) + h*j1/2);
    l2 = X_no(tn(num) + h/2,G_normal(num) + h*k1/2,X_normal(num) + h*l1/2,I_normal(num) + h*j1/2);
    j2 = I_no(tn(num) + h/2,G_normal(num) + h*k1/2,X_normal(num) + h*l1/2,I_normal(num) + h*j1/2,0);

    k3 = G_no(tn(num) + h/2,G_normal(num) + h*k2/2,X_normal(num) + h*l2/2,I_normal(num) + h*j2/2);
    l3 = X_no(tn(num) + h/2,G_normal(num) + h*k2/2,X_normal(num) + h*l2/2,I_normal(num) + h*j2/2);
    j3 = I_no(tn(num) + h/2,G_normal(num) + h*k2/2,X_normal(num) + h*l2/2,I_normal(num) + h*j2/2,0);
    
    k4 = G_no(tn(num) + h,G_normal(num) + h*k3,X_normal(num) + h*l3,I_normal(num) + h*j3);
    l4 = X_no(tn(num) + h,G_normal(num) + h*k3,X_normal(num) + h*l3,I_normal(num) + h*j3);
    j4 = I_no(tn(num) + h,G_normal(num) + h*k3,X_normal(num) + h*l3,I_normal(num) + h*j3,0);

    if num <= x
        G_normal(num+1) = G_normal(num) + h*(k1+2*k2+2*k3+k4)/6;
        X_normal(num+1) = X_normal(num) + h*(l1+2*l2+2*l3+l4)/6;
        I_normal(num+1) = I_normal(num) + h*(j1+2*j2+2*j3+j4)/6;
    
        num = num + 1;
    end
end

plot(tn,G_pat3)
hold on;
plot(tn,G_normal)
xlabel('time')
ylabel("glucose concentration(mg/dl)")
legend({'patient3','normal'},'Location','southwest')
hold on;

u_r = 0:0.0001:2;
g_r = 50:0.0125:300;

mf1_1 = exp(-((u_r-0)/1).^2)
mf1_2 = exp(-((u_r-1)/0.5).^2);
mf1_3 = exp(-((u_r-2)/1).^2);

% figure(2)
% plot(u_r,mf1_1,u_r,mf1_2,u_r,mf1_3);
% title('u(k)-gaussion membership functions')
% xlabel('x')
% legend({'mf1_1','mf1_2','mf1_3'},'Location','southwest')
% hold on;

%...........y(k-1)_mf..............
mf2_1 = exp(-((g_r-50)/100).^2)
mf2_2 = exp(-((g_r-175)/50).^2);
mf2_3 = exp(-((g_r-300)/100).^2);

figure(3)
plot(g_r,mf2_1,g_r,mf2_2,g_r,mf2_3);
title('g(k-1)-gaussion membership functions')
xlabel('g(k-1)')
legend({'mf2_1','mf2_2','mf2_3'},'Location','southwest')
hold on;

%.......patient RFNN........
w_1(1) = 0;
w_2(1) = 0;
w_3(1) = 0;
w_4(1) = 0;
w_5(1) = 0;
w_6(1) = 0;
w_7(1) = 0;
w_8(1) = 0;
w_9(1) = 0;
ETA = 0.5;


num = 0;
for i = tn
    num = num + 1;
    
        o3_1 = min(mf1_1(num),mf2_1(num));
        o3_2 = min(mf1_1(num),mf2_2(num));
        o3_3 = min(mf1_1(num),mf2_3(num));
        o3_4 = min(mf1_2(num),mf2_1(num));
        o3_5 = min(mf1_2(num),mf2_2(num));
        o3_6 = min(mf1_2(num),mf2_3(num));
        o3_7 = min(mf1_3(num),mf2_1(num));
        o3_8 = min(mf1_3(num),mf2_2(num));
        o3_9 = min(mf1_3(num),mf2_3(num));
        
        G_pat3I(num) = (w_1(num)*o3_1 + w_2(num)*o3_2+ w_3(num)*o3_3 + w_4(num)*o3_4 ...
            + w_5(num)*o3_5 + w_6(num)*o3_6 + w_7(num)*o3_7 + w_8(num)*o3_8 + w_9(num)*o3_9)...
            /(o3_1 + o3_2 + o3_3 + o3_4 + o3_5 + o3_6 + o3_7 + o3_8 + o3_9);
        
        w_1(num+1) = w_1(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_1;
        w_2(num+1) = w_2(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_2;
        w_3(num+1) = w_3(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_3;
        w_4(num+1) = w_4(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_4;
        w_5(num+1) = w_5(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_5;
        w_6(num+1) = w_6(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_6;
        w_7(num+1) = w_7(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_7;
        w_8(num+1) = w_8(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_8;
        w_9(num+1) = w_9(num) + ETA * (G_pat3(num)-G_pat3I(num)) * o3_9;

end
%.......normal RFNN........
w_1(1) = 0;
w_2(1) = 0;
w_3(1) = 0;
w_4(1) = 0;
w_5(1) = 0;
w_6(1) = 0;
w_7(1) = 0;
w_8(1) = 0;
w_9(1) = 0;
ETA = 0.5;


num = 0;
for i = tn
    num = num + 1;
    
        o3_1 = min(mf1_1(num),mf2_1(num));
        o3_2 = min(mf1_1(num),mf2_2(num));
        o3_3 = min(mf1_1(num),mf2_3(num));
        o3_4 = min(mf1_2(num),mf2_1(num));
        o3_5 = min(mf1_2(num),mf2_2(num));
        o3_6 = min(mf1_2(num),mf2_3(num));
        o3_7 = min(mf1_3(num),mf2_1(num));
        o3_8 = min(mf1_3(num),mf2_2(num));
        o3_9 = min(mf1_3(num),mf2_3(num));
        
        G_normalI(num) = (w_1(num)*o3_1 + w_2(num)*o3_2+ w_3(num)*o3_3 + w_4(num)*o3_4 ...
            + w_5(num)*o3_5 + w_6(num)*o3_6 + w_7(num)*o3_7 + w_8(num)*o3_8 + w_9(num)*o3_9)...
            /(o3_1 + o3_2 + o3_3 + o3_4 + o3_5 + o3_6 + o3_7 + o3_8 + o3_9);
        
        w_1(num+1) = w_1(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_1;
        w_2(num+1) = w_2(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_2;
        w_3(num+1) = w_3(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_3;
        w_4(num+1) = w_4(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_4;
        w_5(num+1) = w_5(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_5;
        w_6(num+1) = w_6(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_6;
        w_7(num+1) = w_7(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_7;
        w_8(num+1) = w_8(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_8;
        w_9(num+1) = w_9(num) + ETA * (G_normal(num)-G_normalI(num)) * o3_9;

end

figure(4)
plot(tn,G_pat3)
hold on
plot(tn,G_pat3I,'--b');
hold on;
plot(tn,G_normal)
hold on
plot(tn,G_normalI,'--r');
ylabel('glucose concentration')
xlabel('time')




%..........patient3..........
function g_p3 = G_p3(t,G,X,I)
    g_p3 = -0*(G-70)-X*G + 0*exp(-1*t);
end

function x_p3 = X_p3(t,G,X,I)
    x_p3 = -0.02*X + 5.3*10^-6*(I-7);
end

function i_p3 = I_p3(t,G,X,I,u)
    i_p3 = -0.3*(I-7) + u;
end

%..........normal..........
function g_no = G_no(t,G,X,I)
    g_no = -0.0317*(G-70)-X*G + 0*exp(-360*t);
end

function x_no = X_no(t,G,X,I)
    x_no = -0.0123*X + 4.92*10^-6*(I-7);
end

function i_no = I_no(t,G,X,I,u)
    if G-79.0353 > 0
        x = G-79.0353;
    else 
        x = 0;
    end
    i_no = -0.2659*(I-7) + 0.0039*x*t + u;
end


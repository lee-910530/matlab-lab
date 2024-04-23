clc
clear all

%.......data........
h = 0.05;
fin = 500;
tn = 0:h:fin
I_pat1 = zeros(1,fin/h);
X_pat1 = zeros(1,fin/h);
G_pat1 = zeros(1,fin/h);
I_pat1(1) = 50;
X_pat1(1) = 0;
G_pat1(1) = 220; 
I_nor1 = zeros(1,fin/h);
X_nor1 = zeros(1,fin/h);
G_nor1 = zeros(1,fin/h);
I_nor1(1) = 364.8;
X_nor1(1) = 0;
G_nor1(1) = 291.2
%.......runge kutta patient........
x =fin/h;

num = 1;
for i = tn
    k1 = fipat1(tn(num),I_pat1(num));
    k2 = fipat1(tn(num) + h/2,I_pat1(num) + h*k1/2);
    k3 = fipat1(tn(num) + h/2,I_pat1(num) + h*k2/2);
    k4 = fipat1(tn(num) + h,I_pat1(num) + h*k3);
    if num <= x
        I_pat1(num + 1) = I_pat1(num) +h*(k1+2*k2+2*k3+k4)/6;
        num = num + 1;
    end
end

num = 1;
for i = tn
    k1 = fxpat1(tn(num),X_pat1(num),num,I_pat1);
    k2 = fxpat1(tn(num) + h/2,X_pat1(num) + h*k1/2,num,I_pat1);
    k3 = fxpat1(tn(num) + h/2,X_pat1(num) + h*k2/2,num,I_pat1);
    k4 = fxpat1(tn(num) + h,X_pat1(num) + h*k3,num,I_pat1);
    if num <= x
        X_pat1(num + 1) = X_pat1(num) +h*(k1+2*k2+2*k3+k4)/6;
        num = num + 1;
    end
end

num = 1;
for i = tn
    k1 = fgpat1(tn(num),G_pat1(num),num,X_pat1);
    k2 = fgpat1(tn(num) + h/2,G_pat1(num) + h*k1/2,num,X_pat1);
    k3 = fgpat1(tn(num) + h/2,G_pat1(num) + h*k2/2,num,X_pat1);
    k4 = fgpat1(tn(num) + h,G_pat1(num) + h*k3,num,X_pat1);
    if num <= x
        G_pat1(num + 1) = G_pat1(num) +h*(k1+2*k2+2*k3+k4)/6;
        num = num + 1;
    end
end
%.......runge kutta normal........
x =fin/h;

num = 1;
for i = tn
    k1 = finor1(tn(num),I_nor1(num));
    k2 = finor1(tn(num) + h/2,I_nor1(num) + h*k1/2);
    k3 = finor1(tn(num) + h/2,I_nor1(num) + h*k2/2);
    k4 = finor1(tn(num) + h,I_nor1(num) + h*k3);
    if num <= x
        I_nor1(num + 1) = I_nor1(num) +h*(k1+2*k2+2*k3+k4)/6;
        num = num + 1;
    end
end

num = 1;
for i = tn
    k1 = fxnor1(tn(num),X_nor1(num),num,I_nor1);
    k2 = fxnor1(tn(num) + h/2,X_nor1(num) + h*k1/2,num,I_nor1);
    k3 = fxnor1(tn(num) + h/2,X_nor1(num) + h*k2/2,num,I_nor1);
    k4 = fxnor1(tn(num) + h,X_nor1(num) + h*k3,num,I_nor1);
    if num <= x
        X_nor1(num + 1) = X_nor1(num) +h*(k1+2*k2+2*k3+k4)/6;
        num = num + 1;
    end
end

num = 1;
for i = tn
    k1 = fgnor1(tn(num),G_nor1(num),num,X_nor1);
    k2 = fgnor1(tn(num) + h/2,G_nor1(num) + h*k1/2,num,X_nor1);
    k3 = fgnor1(tn(num) + h/2,G_nor1(num) + h*k2/2,num,X_nor1);
    k4 = fgnor1(tn(num) + h,G_nor1(num) + h*k3,num,X_nor1);
    if num <= x
        G_nor1(num + 1) = G_nor1(num) +h*(k1+2*k2+2*k3+k4)/6;
        num = num + 1;
    end
end


figure(1)
yline(70)
hold on
plot(tn,G_pat1,tn,G_nor1);
xlabel('t')
legend({'u(t)','y(t)'},'Location','southwest')
hold on;
%...................................................
u_r = 0:0.00002:2;
g_r = 50:0.0025:300;

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
for i = tn(1:end)
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
        
        G_pat1I(num) = (w_1(num)*o3_1 + w_2(num)*o3_2+ w_3(num)*o3_3 + w_4(num)*o3_4 ...
            + w_5(num)*o3_5 + w_6(num)*o3_6 + w_7(num)*o3_7 + w_8(num)*o3_8 + w_9(num)*o3_9)...
            /(o3_1 + o3_2 + o3_3 + o3_4 + o3_5 + o3_6 + o3_7 + o3_8 + o3_9);
        
        w_1(num+1) = w_1(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_1;
        w_2(num+1) = w_2(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_2;
        w_3(num+1) = w_3(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_3;
        w_4(num+1) = w_4(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_4;
        w_5(num+1) = w_5(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_5;
        w_6(num+1) = w_6(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_6;
        w_7(num+1) = w_7(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_7;
        w_8(num+1) = w_8(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_8;
        w_9(num+1) = w_9(num) + ETA * (G_pat1(num)-G_pat1I(num)) * o3_9;

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
for i = tn(1:end)
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
        
        G_nor1I(num) = (w_1(num)*o3_1 + w_2(num)*o3_2+ w_3(num)*o3_3 + w_4(num)*o3_4 ...
            + w_5(num)*o3_5 + w_6(num)*o3_6 + w_7(num)*o3_7 + w_8(num)*o3_8 + w_9(num)*o3_9)...
            /(o3_1 + o3_2 + o3_3 + o3_4 + o3_5 + o3_6 + o3_7 + o3_8 + o3_9);
        
        w_1(num+1) = w_1(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_1;
        w_2(num+1) = w_2(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_2;
        w_3(num+1) = w_3(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_3;
        w_4(num+1) = w_4(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_4;
        w_5(num+1) = w_5(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_5;
        w_6(num+1) = w_6(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_6;
        w_7(num+1) = w_7(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_7;
        w_8(num+1) = w_8(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_8;
        w_9(num+1) = w_9(num) + ETA * (G_nor1(num)-G_nor1I(num)) * o3_9;

end

figure(4)
plot(tn,G_pat1)
hold on
plot(tn,G_pat1I,'--b');
hold on;
plot(tn,G_nor1)
hold on
plot(tn,G_nor1I,'--r');

%.......patient........
function f_I_pat1 = fipat1(t,I_t);
    f_I_pat1 = -0.3*(I_t-7)+ 0*t;
end

function f_X_pat1 = fxpat1(t,X_t,num,I_pat1);
    f_X_pat1 = -0.02*X_t + 5.3*10^-6*(I_pat1(num)-7);
end

function f_G_pat1 = fgpat1(t,G_t,num,X_pat1);
    f_G_pat1 = -0*(G_t-70)-X_pat1(num)*G_t+1*exp(-100*t);
end
%.......normal........
function f_I_nor1 = finor1(t,I_t);
    f_I_nor1 = -0.2659*(I_t-7)+ 0*t;
end

function f_X_nor1 = fxnor1(t,X_t,num,I_nor1);
    f_X_nor1 = -0.0123*X_t + 4.92*10^-6*(I_nor1(num)-7);
end

function f_G_nor1 = fgnor1(t,G_t,num,X_nor1);
    f_G_nor1 = -0.0317*(G_t-70)-X_nor1(num)*G_t+1*exp(-100*t);
end

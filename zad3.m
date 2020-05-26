clear all;

%%                          Parametry zastosowane w skrypcie

% Dane
K0 = 4.6;
T0 = 5;
T1 = 2.13;
T2 = 4.67;
Tp = 0.5

% Wyznaczone
Kk = 0.46615;
Tk = 20;

% Obliczone
Kr = 0.6*Kk;
Ti = 0.5*Tk;
Td = 0.12*Tk;
Ki = 1/Ti;
Kd = Td;

t = 0:1:500;

% Hdem = tf([1], [1 3 1]);
% Gdem = [1];
% C = pid(Kr, Ti, Td);
% Loop2 = feedback(C*Hdem, Gdem);
% figure(2);
% step(Loop2);
% hold on;
% step(Hdem);
% legend('Odpowiedź układu z pętlą', 'Odpowiedź obiektu');


%%                          Wyznaczanie transmitancji ciągłej

H = tf(K0, [T1*T2 T1+T2 1], 'InputDelay', T0);

%%                          Zdefiniowanie ciągłego regulatora PID

Gc = pid(Kr, Ki, Kd);

%%                          Stworzenie ujemnej pętli sprzężenia zwrotnego

Loop = feedback(Gc*H, [1]);
step(Loop);

%%                          Rysowanie odpowiedzi skokowych

% step(H, t);
% hold on;
% step(Loop, t);
% legend('Odpowiedź skokowa obiektu', 'Odpowiedź skokowa układu');
% hold off;
% print('screeny/zad3zn.png','-dpng','-r400')

%%              Wyznaczanie parametrów regulacji dla regulatora dyskretnego

r2 = (Kk*Td)/Tp;
r1 = -2.95677 %Kk*((Tp)/(2*Ti)-(2)*(Td/Tp)-1);
r0 = 1.62867;

% r0 = Kk*(1 + Tp/(2*Ti) + Td/Tp);
% r1 = Kk*(Tp/(2*Ti) - 2*Td/Tp - 1);
% r2 = Kk*Td/Tp;

% Równanie różnicowe: y(k) = b1(k-1) + b2(k-2) + a1(k-6) + a2(k-7)
b1 = -1.689;
b2 = 0.7105;
a1 = -0.0449;
a2 = -0.04007;

simend = 100;
u(1:7) = 0;
y(1:7) = 0;
yzad(1:9) = 0;
yzad(10:simend) = 1;
e(1:7) = 0;

for k = 8:simend
    y(k) = b1*y(k-1)+b2*y(k-2)+a1*u(k-6)+a2*u(k-7);
    e(k) = yzad(k)-y(k);
    u(k) = r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
end

% stairs(y); 
% hold on;
% stairs(yzad, 'r--');
% xlabel('k');
% legend('odpowiedź skokowa', 'u', 'Location', 'northwest')

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
Ki = (Kr)*1/Ti;
Kd = Kr*Td;

t = 0:1:70;

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

%%                          Rysowanie odpowiedzi skokowych

step(H, t);
hold on;
step(Loop, t);
legend('Odpowiedź skokowa obiektu', 'Odpowiedź skokowa układu');
hold off;
%print('screeny/pidciagly.png','-dpng','-r400')

%%              Wyznaczanie parametrów regulacji dla regulatora dyskretnego

r2 = (Kr*Td)/Tp;
r1 = Kr*((Tp)/(2*Ti)-(2)*(Td/Tp)-1);
r0 = Kr*(1+(Tp/(2*Ti)+(Td/Tp)));


% Równanie różnicowe: y(k) = b1u(k-11) + b0u(k-12) - a1y(k-1) - a0y(k-2)
b1 = 0.05164;
b0 = 0.04608;
a1 = -1.689;
a0 = 0.7105;

simend = 50;
u(1:12) = 0;
y(1:12) = 0;
yzad(1:14) = 0;
yzad(15:simend) = 1;
e(1:12) = 0;

%%                      Symulacja regulatora PID

for k = 13:simend
    y(k) = -a1*y(k-1)-a0*y(k-2)+b1*u(k-11)+b0*u(k-12);
    e(k) = yzad(k)-y(k);
    u(k) = r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
end

 %%                      Rysowanie wykresu funkcji

stairs(y); 
hold on;
stairs(yzad, 'r');
xlabel('k');
legend('Wyjście układu', 'Wyjście zadane', 'Location', 'northwest')
hold off;
%print('screeny/piddyskretny.png', '-dpng', '-r400')

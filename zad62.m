clear all;

%%                          Parametry zastosowane w skrypcie

% Dane
K0nom = 4.6;

K0mnoznik = 1.13

K0 = K0nom * K0mnoznik;
T0nom = 5;
Tmnoznik = [1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2]
T0 = T0nom * Tmnoznik(11);
T1 = 2.13;
T2 = 4.67;
Tp = 0.5;

% Wyznaczone z wykresu
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

% fig = figure;
% step(H, t);
% hold on;
% step(Loop, t);
% legend('Odpowiedź skokowa obiektu', 'Odpowiedź skokowa układu');
% hold off;
%print('screeny/pidciagly.png','-dpng','-r400')

%%              Wyznaczanie parametrów regulacji dla regulatora dyskretnego

r2 = (Kr*Td)/Tp;
r1 = Kr*((Tp)/(2*Ti)-(2)*(Td/Tp)-1);
r0 = Kr*(1+(Tp/(2*Ti)+(Td/Tp)));


% Równanie różnicowe: y(k) = b1u(k-11) + b0u(k-12) - a1y(k-1) - a0y(k-2)
Hd = c2d(H, Tp);  %zdefiniowanie transmitancji dyskretnej dla zad 6
[n d] = tfdata(Hd, 'v'); %ekstrakcja danych z Hd
b1 = n(2);
b0 = n(3);
a1 = d(2);
a0 = d(3);

display(Hd);
x = 22; %najwyzsze opoźnienie
simend = 200;
u(1:x) = 0;
y(1:x) = 0;
yzad(1:x+2) = 0;
yzad(x+3:simend) = 1;
e(1:x) = 0;

%%                      Symulacja regulatora PID

for k = x+1:simend
    y(k) = -a1*y(k-1)-a0*y(k-2)+b1*u(k-(x-1))+b0*u(k-x);
    e(k) = yzad(k)-y(k);
    u(k) = r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
end

 %%                      Rysowanie wykresu funkcji

% fig2 = figure;
% stairs(y); 
% hold on;
% stairs(yzad, 'r');
% xlabel('k');
% legend('Wyjście układu', 'Wyjście zadane', 'Location', 'northwest')
% hold off;
%print('screeny/piddyskretny.png', '-dpng', '-r400')

%%                  Rysowanie krzywej

% Wyznaczone parametry za pomocą kodu powyżej:
Kmnoznik = [1.59 1.52 1.452 1.398 1.35 1.3 1.265 1.22 1.175 1.15 1.12];
% Obliczanie argumentow i wartosci krzywej
[Tcomputed, Kcomputed] = computeCurve(T0, K0, T0nom, K0nom);

fig3 = figure;
plot(Tmnoznik, Kmnoznik);
%axis([0.2 0.4 0 3.5]);
xlabel('T0/T0nom');
ylabel('K0/K0nom');
print('screeny/krzywapid.png', '-dpng', '-r400')
area = computeStabilityArea(K0nom, T0nom, Kmnoznik, Tmnoznik);


function [x, y] = computeCurve(T0, K0, T0nom, K0nom)
    for i = 1:length(T0)
        y(i) = K0(i)/K0nom;
        x(i) = T0(i)/T0nom;
    end
end

function [area] = computeStabilityArea(K0nom, T0nom, Kmnoznik, Tmnoznik)
    area = 0;
    Kmnoznik = Kmnoznik * K0nom;
    Tmnoznik = Tmnoznik * T0nom;
    for i = 1:length(Kmnoznik)-1
        %Triangle:
       heightT =  Kmnoznik(i) - Kmnoznik(i+1);
       widthT = Tmnoznik(i+1) - Tmnoznik(i);
       triangleArea = widthT*heightT*0.5;
       %Square
       heightS = Kmnoznik(i+1);
       widthS = widthT;
       squareArea = heightS * widthS;
       area = area + triangleArea + squareArea;
    end
end







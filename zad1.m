%%                            Wyznaczanie danych



%% Parametry zastosowane w skrypcie
K0 = 4.6;
T0 = 5;
T1 = 2.13;
T2 = 4.67;
Tp = 0.5;

%% Wyznaczenie licznika i mianownika dla funkcji c2d
H = tf([4.6],[T1*T2 T1+T2 1], 'InputDelay', T0); %transmitancja ciągła

%% Obliczenie licznika i mianownika transmitancji dyskretnej

% Funkcja c2d pozwala na wyznaczenie transmitancji dyskretnej z wybranym
% ekstrapolatorem (domyślnie zoh)
[N, D] = c2d(H, Tp); % ekstrapolator zerowego rzędu
Hd = tf(N, D, Tp); %transmitancja dyskretna

%%                           Rysowanie wykresów:
hold on;
figure(1);




%% Wyznaczenie odpowiedzi skokowej transmitancji ciągłej
step(H);
stepI = stepinfo(H); %parametry tej odpowiedzi

%% Wyznaczenie odpowiedzi skokowej transmitancji dyskretnej
step(Hd);
stepII = stepinfo(Hd);




grid on;
xlabel("t");
ylabel("y(t)");
title("Porównanie odpowiedzi skokowych obu transmitancji");
hold off;
print('odps.png','-dpng','-r400')

%%                   Porównanie współczynników wzmocnienia statycznego

%% Transmitancja ciągła

s = 0;
Hstat = exp(-5*s)*(4.6/(9.947*s^2+6.8*s+1));

%% Transmitancja dyskretna

z = 1;
Hdstat = (0.05164*z + 0.04608)/(z^2 - 1.689*z + 0.7105);
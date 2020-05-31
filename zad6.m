clear all;

%%                  Definicja parametrow uzywanych w skrypcie

% DMC
% Wspolczynniki transmitancji dyskretnej
a1 = -1.689;
a0 = 0.7105;
b1 = 0.05164;
b0 = 0.04608;

% Stale wartosci
Lm = [b1 b0];
Mm = [1 a1 a0 0 0 0 0 0 0 0 0 0 0];
simtime = 130; 
tfunc = tf(Lm,Mm,0.5);
timebase = 0:0.5:simtime
s = step(tfunc,timebase);
changetime = 15;
J = 0;

% Parametry do dopasowania
D = 95;
N = 18;
Nu = 6;
lambda = 4;
step = 1;

% PID
% Dane
K0 = 4.6;
T0 = 5;
T1 = 2.13;
T2 = 4.67;
Tp = 0.5;
Jpid = 0;

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

%%              Przygotowanie macierzy oraz sygnalow

u(1:12) = 0;
ymod(1:12) = 0;
yzad(1:14) = 0;
yzad(changetime:simtime) = 1;
dupk = zeros(1, D-1);
M=zeros(N,Nu);
Mp=zeros(N,D-1);

M = fillM(N, Nu, M, s);
Mp = fillMp(N, D, Mp, s);

%%                       Wyznaczanie macierzy K

I = eye(Nu);
K = (M'*M+lambda*I)^-1*M';
Ku = K(1,:)*Mp;
Ke = sum(K(1,:));

%%                      Symulacja algorytmu DMC

for k = 13:simtime
    ymod(k) = -a1*ymod(k-1)-a0*ymod(k-2)+b1*u(k-11)+b0*u(k-12);
    for j = D-1:-1:2
        dupk(j) = dupk(j-1);
    end
    e(k) = yzad(k)-ymod(k);
    descriptionForY = sprintf("k=%d,e(k)=%0.5f",k,e(k));
    duk = Ke*e(k)-Ku*dupk';
    dupk(1) = duk; 
    u(k) = u(k-1)+duk; 
    J = quality(J, k, yzad, ymod);
    
end

regulationTimeDMC = stablePoint(e) - changetime;

%%              Wyznaczanie parametrów regulacji dla dyskretnego PID

r2 = (Kr*Td)/Tp;
r1 = Kr*((Tp)/(2*Ti)-(2)*(Td/Tp)-1);
r0 = Kr*(1+(Tp/(2*Ti)+(Td/Tp)));

% Równanie różnicowe: y(k) = b1u(k-11) + b0u(k-12) - a1y(k-1) - a0y(k-2)
b1 = 0.05164;
b0 = 0.04608;
a1 = -1.689;
a0 = 0.7105;

u2(1:12) = 0;
y(1:12) = 0;
clear yzad;
yzad(1:14) = 0;
yzad(15:simtime) = 1;
clear e;
e(1:12) = 0;

%%                      Symulacja regulatora PID

for k = 13:simtime
    y(k) = -a1*y(k-1)-a0*y(k-2)+b1*u2(k-11)+b0*u2(k-12);
    e(k) = yzad(k)-y(k);
    u2(k) = r2*e(k-2)+r1*e(k-1)+r0*e(k)+u2(k-1);
    Jpid = quality(J, k, yzad, y);
end

regulationTimePID = stablePoint(e) - changetime;

%%                                  Rysowanie wykresow

descriptionForY = sprintf("Porównanie działania obu metod regulacji");
descriptionForU = sprintf("Porównanie sygnałów sterujących");
fig = figure;
time = 0:step:simtime-1;

% Y
subplot(2, 1, 1);
stairs(time, ymod);
hold on;
stairs(time, y);
stairs(time, yzad, 'k--');
preparePlot(simtime, descriptionForY, 1.5);
leg = legend('DMC', 'PID', 'Yzad');
set(leg, 'location', 'best');
hold off;

% U
subplot(2, 1, 2);
stairs(time, u);
hold on;
stairs(time, u2);
preparePlot(simtime, descriptionForU, 2.0);
legend('DMC', 'PID');
hold off;

print('screeny/zad6.png','-dpng','-r400');

%%                                              Funkcje

function [] = preparePlot(simtime, description, yAxHeight)
    axis([0 simtime 0 yAxHeight]);
    xlabel('k');
    ylabel('y(k), yzad(k)');
    title(description);
    grid on;
end


function [point] = stablePoint(s)
    differenceReq = 0.01;
    lenReq = 5;
    result = 0;
    for i = 26:length(s)-lenReq
        for j = 0:lenReq-1
            x = s(i+j)
            if abs(x) < differenceReq
                result = result + 1;
            end
        end
        if result == lenReq
            point = i;
            break;
        end
        result = 0; 
    end
end

function [J] = quality(J, k, yzad, ymod)
    J = J+(yzad(k)-ymod(k))^2;
end


function [M] = fillM(N, Nu, M, s)
    for i = 1:N
        for j = 1:Nu
            if(i >= j)
                M(i,j) = s(i-j+1);
            end
        end
    end
end

function [Mp] = fillMp(N, D, Mp, s)
    for i = 1:N
        for j = 1:D-1
            if i+j <= D
                Mp(i,j) = s(i+j)-s(j);
            else
                Mp(i,j) = s(D)-s(j);
            end
        end
    end
end


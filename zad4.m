clear all;

%%                  Definicja parametrow uzywanych w skrypcie

% Wspolczynniki transmitancji dyskretnej
a1 = -1.689;
a0 = 0.7105;
b1 = 0.05164;
b0 = 0.04608;

% Stale wartosci
Lm = [b1 b0];
Mm = [1 a1 a0 0 0 0 0 0 0 0 0 0 0];
simtime = 55; 
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
    opis = sprintf("k=%d,e(k)=%0.5f",k,e(k));
    duk = Ke*e(k)-Ku*dupk';
    dupk(1) = duk; 
    u(k) = u(k-1)+duk; 
    J = quality(J, k, yzad, ymod);
    
end

%%                      Sprawdzenie jakosci regulacji i czasu regulacji

regulationTime = stablePoint(e) - changetime;
display(J);
display(regulationTime);
display(regulationTime+changetime);

%%                                  Rysowanie wykresow

str = sprintf("D=%d,N=%d,Nu=%d,lambda=%0.1f: J=%d t=%d",D,N,Nu,lambda, J, regulationTime)
fig = figure;
czas = 0:step:simtime-1;      
subplot(2,1,1);
stairs(czas, ymod)
axis([0 simtime 0 1.5])
xlabel('k');
title(str)
hold on;
grid on;
ylabel('y(k),yzad(k)')
stairs(1:simtime,yzad,'k--')

subplot(2,1,2);
stairs(czas, u)
grid on;
xlabel('k');
ylabel("u(k)")
title("Sygnal sterujacy u(k)")
set(gcf,'Units','centimeters','Position', [0 0 14 10]);
print('screeny/podsumowanie.png','-dpng','-r400');

%%



%%                                              Funkcje

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





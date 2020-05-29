J = [14.3 14.42 14.56 14.63 14.64 14.64];
t = [37 36 31 31 31 31];
wage = 30;

scatter(J, t);
w = checkBest(J, t, wage);
print('screeny/5bscatter.png','-dpng','-r400');
display(w)

function [w] = checkBest(J, t, wage)
w = [0 0];
c = 10000;
for i = 1:length(J)
    squaresSum = wage*J(i)^2 + t(i)^2;
    module = sqrt(squaresSum);
    if module < c
        c = module;
        w(1) = J(i);
        w(2) = t(i);
    end
end
end


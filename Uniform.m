clear;
close all;

% Generator liczb (pseudo)losowych oparty na przekształceniu piłokształtnym
function PRN = rand_p(rows, cols, z)
    persistent X_pom;
    if isempty(X_pom)
       X_pom = cputime * 0.736217631286;
    end
    X_0 = X_pom;
    
    PRN = zeros(rows, cols);
    for i=1:1:rows
       for j=1:1:cols
            PRN(i, j) = X_0 * z - floor(X_0 * z);
            X_0 = PRN(i, j);
       end
    end
    X_pom = X_pom*1.01;
end

% Liniowy generator liczb (pseudo)losowych - LCG
function PRN = rand_l(rows, cols, m, c)
    persistent X_pom;
    if isempty(X_pom)
       X_pom = cputime * 0.736217631286;
    end
    
    a_0 = 1997;
    PRN = zeros(rows, cols);
    A = [];
    
    for k=1:1:(rows*cols)
        A(end+1) = a_0;
        a_0 = a_0 + 2;
    end
    
    t = X_pom:X_pom:rows*cols*X_pom;   
    X = cos(2*pi*1000*t);
    
    for i=1:1:rows
       for j=1:1:cols
           PRN(i, j) = mod(dot(A, X) + c, m);
           X(end) = [];
           X = [PRN(i, j),X];
       end
    end
    X_pom = X_pom*1.01;
end

% Badanie okresu powtarzania się wygenerowanych próbek
function Int = CheckInterval(X)
    for i=1:1:length(X)
        Int = 0;
       for j=i:1:length(X)
          if((X(i) == X(j)) && (i~=j))
              return;
          end
          Int = Int + 1;
       end
    end
    Int = 0;
end

% Testowa metoda Monte Carlo wyznaczająca wartość pi
function pi_est = MonteCarlo(X, Y)
    
    pkt_w_kwadracie = length(X);
    pkt_w_kole = 0;
    
    for i=1:1:length(X)
       if(X(i)*X(i) + Y(i)*Y(i) <= 1)
          pkt_w_kole = pkt_w_kole + 1; 
       end
    end
    pi_est = 4 * (pkt_w_kole / pkt_w_kwadracie);
end

% Metoda do sprawozdania
function pie = DrawMonteCarlo(X, Y)
    x = 0;
    y = 0;
    r = 1;
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    
%     X = rand_p(1,1000,5);
%     Y = rand_p(1,1000,5);
    
    X = 2*X-1;
    Y = 2*Y-1;
    
    pie = MonteCarlo(X, Y);
    figure(1);
    plot(x+xp,y+yp);
    hold on;
    plot(X, Y, '.');
    title('Metoda Monte Carlo, przyblizenie wielkosci liczby pi','interpreter','latex');
    legend('Okrąg jednostkowy','Wygenerowane próbki');
end

% Metoda do sprawozdania
function DrawHistogram()
    X1 = rand_p(1,10000,5);
    X2 = rand_l(1,10000,1,1);
    X3 = rand(1,10000);
    
    subplot(3,1,1);
    histogram(X1);
    title('Generator liczb (pseudo)losowych oparty na przeksztalceniu piloksztaltnym - N=10000','interpreter','latex')
    subplot(3,1,2);
    histogram(X2);
    title('Liniowy generator liczb (pseudo)losowych - LCG - N=10000','interpreter','latex');
    subplot(3,1,3);
    histogram(X3);
    title('Generator liczb (pseudo)losowych rand - wbudowany w oprogramowanie Matlab - N=10000','interpreter','latex');
end

% Metoda do sprawozdania
function Complexity()
    N = 1000:1000:20000;
    T = zeros(1,length(N));
    for i=1:1:length(N)
        tic
        rand(1, N(i));
        T(i) = toc;
    end
    figure(1);
    plot(N, T, '.');
    title('Zlozonosc czasowa','interpreter','latex');
    xlabel('Ilosc probek - N','interpreter','latex');
    ylabel('Czas wykonania algorytmu - t[s]','interpreter','latex');
end

% Test par
function [h, prc, avg] = ChiSqrTest(W, m)
    X = zeros(1,length(W)/2);
    Y = zeros(1,length(W)/2);
    pom = 1;
    for t=1:2:length(W)-1
        X(pom) = W(t);
        pom = pom + 1;
    end
    pom = 1;
    for z=2:2:length(W)
        Y(pom) = W(z);
        pom = pom + 1;
    end
    h = zeros(m, m);
    for i=1:1:m
       for j=1:1:m
           for k=1:1:length(X)
              if(i == 1 && j == 1)
                 if(X(k) <= i/m && Y(k) <= j/m)
                    h(i, j) = h(i, j) + 1;
                 end
              elseif(i == 1 && j ~= 1)
                  if((X(k) <= i/m && Y(k) <= j/m) && (Y(k) > (j-1)/m))
                    h(i, j) = h(i, j) + 1;
                  end
              elseif(i ~= 1 && j == 1)
                  if((X(k) <= i/m && Y(k) <= j/m) && (X(k) > (i-1)/m))
                    h(i, j) = h(i, j) + 1;
                  end
              elseif(i ~= 1 && j ~= 1)
                  if((X(k) <= i/m && Y(k) <= j/m) && (X(k) > (i-1)/m && Y(k) > (j-1)/m))
                    h(i, j) = h(i, j) + 1; 
                  end
              end
           end
       end
    end
    prfct = (length(W)/2)/m^2;
    prc = zeros(m,m);
    for l=1:1:m
       for n=1:1:m
           prc(l, n) = abs(100 - (h(l, n)/prfct * 100));
       end
    end
    avg = sum(sum(prc))/m^2;
end

% Gra w chaos
function Sierpinski()
    Xp = [2, -2, 0];
Yp = [0, 0, 4];
Ax = 2;
Ay = 0;
Bx = -2;
By = 0;
Cx = 0;
Cy = 4;
pX = 0;
pY = 2;

for i=1:1:100000
    Xp(end+1) = pX;
    Yp(end+1) = pY;
    a = rand_p(1,1,5);
    if(a > 0 && a <= 0.33333)
       pX = pX + (Ax - pX) / 2;
       pY = pY + (Ay - pY) / 2;
    elseif(a > 0.33333 && a <= 0.66666)
       pX = pX + (Bx - pX) / 2;
       pY = pY + (By - pY) / 2;
    elseif(a > 0.66666 && a <= 1)
       pX = pX + (Cx - pX) / 2;
       pY = pY + (Cy - pY) / 2;
    end
end
plot(Xp, Yp, '.');
title('Gra w chaos - funkcja rand','interpreter','latex');
xlabel('X','interpreter','latex');
ylabel('Y','interpreter','latex');
end

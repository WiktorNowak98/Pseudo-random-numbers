clear;
close all;

% x = -80:0.01:80;
% y = 0.5*exp(-abs(x));
% % x = 0:0.001:1;
% % y = 2*x;
% 
% new_Y = calkowanie(x, y);
% Rozklad = odwr(x, new_Y);
% 
% y_moj = dens_func(x, Rozklad);
% 
% rozklad = Zad_4();
% y_moj2 = dens_func(x, rozklad); 
% 
% % histogram(Rozklad);
% % title('Numeryczny algorytm wyznaczenia rozkladu - rozklad z zadania 1','interpreter','latex');
% 
% Blad_numerycznego = MSE(y_moj, y);
% Blad_analitycznego = MSE(y_moj2, y);
% 
% figure(1);
% plot(x,y);
% hold on;
% plot(x,y_moj);
% hold on;
% plot(x, y_moj2);
% legend('Gestosc rozkladu','Gestosc numeryczna','Gestosc analityczna', 'interpreter', 'latex');
% title('Porownanie jakosci rozkladow uzyskanych roznymi metodami', 'interpreter', 'latex');
% xlabel('X','interpreter', 'latex');
% ylabel('Y','interpreter','latex');

[xp, xl, xr] = Zad_4();

histogram(xp);
%     subplot(3,1,1);
%     histogram(xp);
%     title('Generator oparty na generatorze piloksztaltnym i danym rozkladzie','interpreter','latex')
%     subplot(3,1,2);
%     histogram(xl);
%     title('Generator oparty na generatorze liniowym i danym rozkladzie','interpreter','latex');
%     subplot(3,1,3);
%     histogram(xr);
%     title('Generator oparty na funkcji rand i danym rozkladzie','interpreter','latex');
% 
% iteracje = 10;
% x = 0:0.01:1;
% y = 2*x;
% 
% [rp, rl, rr] = Zad_1();
% yp = dens_func(x,rp);
% yl = dens_func(x,rl);
% yr = dens_func(x,rr);
% 
% figure(1);
% plot(x,y);
% hold on;
% plot(x,yp);
% hold on;
% plot(x,yl);
% hold on;
% plot(x, yr);
% legend('Oryginal','Piloksztaltny','Liniowy','Rand','interpreter','latex');
% title('Porownanie oryginalnej funkcji gestosci z przyblizonymi za pomoca generatorow','interpreter','latex');
% xlabel('X','interpreter','latex');
% ylabel('Y','interpreter','latex');

% Błąd średniokwadratowy MSE
function x = MSE(y2, y3)
    suma = 0;
    for i=1:1:length(y2)
        a = (y2(i) - y3(i))^2;
        suma  = suma + a;
    end
    x = suma/length(y2);
end

% Gestosc
function y = dens_func(X, Y)
    pd = fitdist(Y', 'Kernel', 'Kernel', 'epanechnikov');
    y = pdf(pd, X);
end

% Calkowanie numeryczne - suma Riemanna
function new_Y = calkowanie(X, Y)
    new_Y = [];
    x0 = 0;
    suma = 0;
    for i=1:1:length(X)
        suma = Y(i)*(abs(X(i) - x0)) + suma;
        new_Y(end+1) = suma;
        x0 = X(i);
    end
end

% Numeryczne odwracanie dystrybuanty
function Rozklad = odwr(X, Y)
    new_X = Y;
    new_Y = X;
    iterations = 100000;
    
    los_X = rand_l(1,iterations,1,1);
    Pom = [];
    Rozklad = [];
    
    start = 1000000;
    gdzie = 0;
    for i=1:1:iterations
        gdzie = 0;
        start  = 10000000;
       for j=1:1:length(new_X)
           if(abs(los_X(i) - new_X(j)) < start)
               start = abs(los_X(i) - new_X(j));
               gdzie = j;
           end  
       end
       Rozklad(end+1) = new_Y(gdzie);
    end
end

% Rozklad |2x
function [Rozkladp, Rozkladl, Rozkladr] = Zad_1()
    X1 = rand_p(1,100000,5);
    X2 = rand_l(1,100000,1,1);
    X3 = rand(1,100000);
    Rozkladp = [];
    Rozkladl = [];
    Rozkladr = [];
    
    for i=1:1:length(X1)
        Rozkladp(end+1) = sqrt(X1(i));
    end
    
    for i=1:1:length(X1)
        Rozkladl(end+1) = sqrt(X2(i));
    end
    
    for i=1:1:length(X1)
        Rozkladr(end+1) = sqrt(X3(i));
    end
    
end

% Rozklad |x+1, -x+1
function [Rozkladp, Rozkladl, Rozkladr] = Zad_2()
    X1 = rand_p(1,100000,5);
    X2 = rand_l(1,100000,1,1);
    X3 = rand(1,100000);
    Rozkladp = [];
    Rozkladl = [];
    Rozkladr = [];
    
    for i=1:1:length(X1)
       if(X1(i) >= 0 && X1(i) < 1/2)
           Rozkladp(end+1) = sqrt(2*X1(i))-1;
       elseif(X1(i) >= 1/2 && X1(i) < 1)
           Rozkladp(end+1) = 1 - sqrt(2-2*X1(i));
       end
    end
    
    for i=1:1:length(X2)
       if(X2(i) >= 0 && X2(i) < 1/2)
           Rozkladl(end+1) = sqrt(2*X2(i))-1;
       elseif(X2(i) >= 1/2 && X2(i) < 1)
           Rozkladl(end+1) = 1 - sqrt(2-2*X2(i));
       end
    end
    
    for i=1:1:length(X3)
       if(X3(i) >= 0 && X3(i) < 1/2)
           Rozkladr(end+1) = sqrt(2*X3(i))-1;
       elseif(X3(i) >= 1/2 && X3(i) < 1)
           Rozkladr(end+1) = 1 - sqrt(2-2*X3(i));
       end
    end
end

% Rozklad wykładniczy
function [Rozkladp, Rozkladl, Rozkladr] = Zad_3()
    X1 = rand_p(1,100000,5);
    X2 = rand_l(1,100000,1,1);
    X3 = rand(1,100000);
    Rozkladp = [];
    Rozkladl = [];
    Rozkladr = [];
    
    for i=1:1:length(X1)
       Rozkladp(end+1) = -log(1-X1(i)); 
    end
    
    for i=1:1:length(X2)
       Rozkladl(end+1) = -log(1-X2(i)); 
    end
    
    for i=1:1:length(X3)
       Rozkladr(end+1) = -log(1-X3(i)); 
    end
end

% Rozkład Laplace'a (oszukany)
function [Rozkladp, Rozkladl, Rozkladr] = Zad_4()
    X1 = rand_p(1,100000,5);
    X2 = rand_l(1,100000,1,1);
    X3 = rand(1,100000);
    Rozkladp = [];
    Rozkladl = [];
    Rozkladr = [];
    z = 0;
    
    for i=1:1:length(X1)
        los = rand(1,1);
        if(los >= 0 && los < 1/2)
            z = -1;
        elseif(los >= 1/2 && los < 1)
            z = 1;
        end
        Rozkladp(end+1) = -log(1-X1(i))*z;
    end
    
    z = 0;
    
    for i=1:1:length(X2)
        los = rand(1,1);
        if(los >= 0 && los < 1/2)
            z = -1;
        elseif(los >= 1/2 && los < 1)
            z = 1;
        end
        Rozkladl(end+1) = -log(1-X2(i))*z;
    end
    z = 0;
    
    for i=1:1:length(X3)
        los = rand(1,1);
        if(los >= 0 && los < 1/2)
            z = -1;
        elseif(los >= 1/2 && los < 1)
            z = 1;
        end
        Rozkladr(end+1) = -log(1-X3(i))*z;
    end
end

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
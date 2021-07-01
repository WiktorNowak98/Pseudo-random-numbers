clear;
close all;

f = @(x) ((-1<x) && (x<=0)) .*(x+1)+((0<x) && (x<1)).*(-x+1);
Rozklad = Generuj_liczby(f, 1000, -1, 1);
Rozklad1 = Generuj_liczby(f,100000,-1,1);
Rozklad2 = Generuj_liczby(f,100000,-0.5,0);
Rozklad3 = Generuj_liczby(f,100000,-0.25,0.75);

figure(1);
subplot(4,1,1);
histogram(Rozklad,'Normalization','pdf');
title('Generator liczb losowych z dowolnego rozkladu - Zadanie 1 n = 1000 - [a,b] = [-1,1]','interpreter','latex');
subplot(4,1,2);
histogram(Rozklad1,'Normalization','pdf');
title('Generator liczb losowych z dowolnego rozkladu - Zadanie 1 n = 100000 - [a,b] = [-1,1]','interpreter','latex');
subplot(4,1,3);
histogram(Rozklad2,'Normalization','pdf');
title('Generator liczb losowych z dowolnego rozkladu - Zadanie 1 n = 100000 - [a,b] = [-0.5,0]','interpreter','latex');
subplot(4,1,4);
histogram(Rozklad3,'Normalization','pdf');
title('Generator liczb losowych z dowolnego rozkladu - Zadanie 1 n = 100000 - [a,b] = [-0.25,0.75]','interpreter','latex');

function y = dens_func(X, Y)
    pd = fitdist(Y', 'Kernel', 'Kernel', 'epanechnikov');
    y = pdf(pd, X);
end

function Rozklad = Zad_1()
    
    fun = @(x)((-1<x) && (x<=0)) .*(x+1)+((0<x) && (x<1)).*(-x+1);
    U1 = 2*(rand(1,100000))-1;
    U2 = 1*rand(1,100000);
    
    Rozklad = [];
    for i=1:1:length(U1)
        if(U2(i) <= fun(U1(i)))
            Rozklad(end+1) = U1(i);
        end
    end
end

function Rozklad = Zad_2()
    fun = @(x) ((x>0) && (x<=1/100)) .* (50) + ((x>1/100) && (x<=1)) .* (100/198);
    U1 = rand(1,100000);
    U2 = 50*rand(1,100000);
    
    Rozklad = [];
    for i=1:1:length(U1)
       if(U2(i) <= fun(U1(i)))
          Rozklad(end+1) = U1(i); 
       end
    end
end

function Rozklad = Zad_2_lepiej()
    fun = @(x) ((x>0) && (x<=1/100)) .* (50) + ((x>1/100) && (x<=1)) .* (100/198); % Szukamy
    gun = @(x) ((x>0) && x<1).*(-2*x+2);
    c = 25;
    
    X1 = rand(1,100000);
    V = [];
    
    for i=1:1:length(X1)
        V(end+1) = 1 - sqrt(1-X1(i));
    end
    U = rand(1,length(V));
    Rozklad = [];
    for i=1:1:length(U)
       if(U(i) * c * gun(V(i)) <= fun(V(i)))
           Rozklad(end+1) = V(i);
       end
    end
end

function Rozklad = Zad_2_exp()
    fun = @(x) ((x>0) && (x<=1/100)) .* (50) + ((x>1/100) && (x<=1)) .* (100/198); % Szukamy
    gun = @(x) exp(-x);
    X1 = rand(1,100000);
    V = [];
    c = 50;
    
    for i=1:1:length(X1)
       V(end+1) = -log(1-X1(i));  
    end
    %histogram(V);
    U = rand(1,length(V));
    Rozklad = [];
    for i=1:1:length(U)
       if(U(i)* c * gun(V(i)) <=  fun(V(i)))
          Rozklad(end+1) = V(i); 
       end
    end
end

function Rozklad = Zad_3()
    r = sqrt(2/pi);
    fun = @(x) ((x>=-r) && (x<=r)) .* (sqrt(r^2 - x^2));
    
    U1 = 2*r*rand(1,100000)- r; % od -r do r
    U2 = rand(1,100000)*r; % Y od 0 do 1
    
    Rozklad = [];
    for i=1:1:length(U1)
       if(U2(i) <= fun(U1(i)))
           Rozklad(end+1) = U1(i);
       end
    end
end

function Rozklad = Zad_4()
    fun = @(x) (1/sqrt(2*pi))*exp(-x^2/2);
    gun = @(x) 0.5*exp(-abs(x));
    
    c = sqrt(2*2.71828182846/pi);
    U = rand(1,100000);
    X3 = rand(1,100000);
    V = [];
    z = 0;
    
    for i=1:1:length(X3)
        los = rand(1,1);
        if(los >= 0 && los < 1/2)
            z = -1;
        elseif(los >= 1/2 && los < 1)
            z = 1;
        end
        V(end+1) = -log(1-X3(i))*z;
    end
    Rozklad = [];
    for i=1:1:length(U)
       if(U(i) * c * gun(V(i)) <= fun(V(i)))
           Rozklad(end+1) = V(i);
       end
    end
end

function Rozklad = Box_Muller(iter, mu, sigma)
    
    Rozklad = [];
    for i=1:1:iter/2
        u1 = rand(1,1);
        u2 = rand(1,1);
        Pierw = sigma * sqrt(-2*log(u1));
        Rozklad(end+1) = Pierw * cos(2*pi*u2) + mu;
        Rozklad(end+1) = Pierw * sin(2*pi*u2) + mu;
    end
end

function Rozklad = Generuj_liczby(fun, iter, a, b)
    x = a:0.01:b;
    y = zeros(1,length(x));
    for i=1:1:length(x)
        y(i) = fun(x(i));
    end
    d = ceil(max(y));
    Rozklad = [];
    
    going = true;
    while(going)
        U1 = a + (b - a) * rand(1,1);
        U2 = d * rand(1,1);
        if(U2 <= fun(U1))
            Rozklad(end+1) = U1;
        end
        if(length(Rozklad) == iter)
            going = false;
        end
    end
end

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

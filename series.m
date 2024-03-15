%% Cargar modulos solo Octave
% pkg load symbolic
clear
clc
close all
syms x real
% Progama que calcula las series de Fourier con n terminos

% Ingresar datos
L = sym(pi)                                             % Los intervalos van de [-L a L]
n = 10;

%f(x) = x^2;                                            % Definir la funcion periódica
%f(x) = piecewise(x<0, -1, x>0, 1)                       % Función escalón en dos intervalos
f(x) = piecewise(x<0, 0, x>0, x)
f_x = f(x);
f_minus_x = f(-x);

% Definir las variables para el calculo de las series
a0 = (1/(2*L)) * int(f, -L, L);
an = sym(0).*ones(1, n);
bn = sym(0).*ones(1, n);
Fourier_series = a0;

% Determinar si la funcion es par o impar y calculo de la serie de Fourie
% Calculo de las series y Grafica
figure;
hold on;

if simplify(f_x - f_minus_x) == 0
  disp("La función es par")
  for k = 1:n
    an(k) = (1/L) * int(f * cos((k * sym(pi) * x) / L), -L, L);
    Fourier_series = Fourier_series + an(k) * cos((k * sym(pi) * x) / L);
    serie = a0 + an(k)*cos((k * sym(pi) * x) / L);
    % ezplot(serie,[-L,L])
    plot_handle = ezplot(serie);
    set(plot_handle, "Color", [65/255, 43/255, 21/255]);
  end


elseif simplify(f_x + f_minus_x) == 0
  disp("La función es impar")
  for k = 1:n
    bn(k) = (1/L) * int(f * sin((k * sym(pi) * x) / L), -L, L);
    Fourier_series = Fourier_series + bn(k) * sin((k * sym(pi) * x) / L);
    serie = a0 + bn(k)*sin((k * sym(pi) * x) / L);
    % ezplot(serie,[-L,L],[1 0 0])
    plot_handle = ezplot(serie);
    set(plot_handle, "Color", [0 1 1]);
  end

else
  disp("La función no es par ni impar")
  for k = 1:n
    an(k) = (1/L) * int(f * cos((k * sym(pi) * x) / L), -L, L);
    bn(k) = (1/L) * int(f * sin((k * sym(pi) * x) / L), -L, L);
    serie_an = a0 + an(k) * cos((k * sym(pi) * x) / L)
    serie_bn = a0 + bn(k) * sin((k * sym(pi) * x) / L)
    Fourier_series = Fourier_series + an(k) * cos((k * sym(pi) * x) / L) + bn(k) * sin((k * sym(pi) * x) / L);
    plot_handle = ezplot(serie_an);
    set(plot_handle, "Color", [65/255, 43/255, 21/255]);
    plot_handle = ezplot(serie_bn);
    set(plot_handle, "Color", [0 1 1]);
  end
end

x = -double(L):0.1:double(L);
y = double(f(x));

%ezplot(Fourier_series)
plot_handle = ezplot(Fourier_series);
set(plot_handle, "Color", [0 0 1], "linewidth", 1.5);
plot(x,y,'red','Linewidth',1.5)
xlim([-double(L), double(L)])
grid on
hold off
title('Series de Fourier Analisis');
legend('Coeficiente a_n','Coficiente b_n');
xlabel('x');
ylabel('f(x)');

# Grafica 2
figure
hold on;
% ezplot(Fourier_series)
plot_handle = ezplot(Fourier_series);
set(plot_handle, "Color", [0 0 1], "linewidth", 1.5);
plot(x,y,'red','Linewidth',1.5)
xlim([-double(L), double(L)])
grid on
hold off
title('Series de Fourier Resultado');
legend('Funcion original','Serie de fourier');
xlabel('x');
ylabel('f(x)');

% Cálculo de los coeficientes de la serie de Fourier
%{
for k = 1:n
    an(k) = (1/L) * int(f * cos((k * sym(pi) * x) / L), -L, L);
    bn(k) = (1/L) * int(f * sin((k * sym(pi) * x) / L), -L, L);
end
%}

% Creación de la función de la serie de Fourier
%{
for k = 1:n
    Fourier_series = Fourier_series + an(k) * cos((k * sym(pi) * x) / L) + bn(k) * sin((k * sym(pi) * x) / L);
end
%}

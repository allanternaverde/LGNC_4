clc
clear all
close all

%%Pouso lunar
%Dados do problema
g = -1.5;
E = 2.5;
theta = 35*pi/180;
m = 1;

%% Plotando as curvas y - y', para altitude m�xima 500 m
% Descobrindo as v. iniciais para pauso suave com uma config. de propulsor 
v_0_1 = -sqrt(2*((E/m) + g)*(500));
v_0_2 = -sqrt(2*((E*cos(theta)/m) + g)*(500));

% Plotando os gr�ficos de pouso suave
for i = 1:500001
    y(i) = 500 - (0.001*(i-1));
    y_linha1(i) = -sqrt(((v_0_1)^2) + (2*((E/m) + g)*(y(i) - 500)));
    y_linha2(i) = -sqrt(((v_0_2)^2) + (2*((E*cos(theta)/m) + g)*(y(i) - 500)));
end

%Plotando o ponto de condi��o inicial do problema
y_inicial = 150;
y_linha_inicial = -14;
x_inicial = 0;
x_linha_inicial = 0;


figure(1)
plot(y_linha1,y,y_linha2,y,y_linha_inicial,y_inicial,'*');
grid
title('Regi�o de Pouso Suave')
ylabel('Altitude (m)')
xlabel('Velocidade (m/s)')

clear y
clear y_linha1

%% Definindo pouso suave em Y
% Estrat�gia 1: queda livre + propulsor 0�.

%Altitude para ligar o propulsor
y_1 = ((y_linha_inicial^2) - (2*g*y_inicial))/(2*(E/m));
fprintf('Estrat�gia 1 - Queda Livre + Propulsor a 0�: \n')
fprintf('A altitude de acionamento do propulsor � %.2f m.\n\n',y_1)
passo = floor((y_inicial-y_1)/0.001);
passo_1 = floor(y_1/0.001);

for i = 1:(passo+1)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha_inicial)^2) + (2*g*(y(i) - y_inicial)));
end

for i = (passo+1):(passo + passo_1 + 2)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha1(passo + 1))^2) + (2*((E/m) + g)*(y(i) - y(passo + 1))));
end

figure(2)
plot(y_linha1,y)
grid
title('Estrat�gia 1 - Queda Livre + Propulsor a 0�')
ylabel('Altitude (m)')
xlabel('Velocidade (m/s)')


clear y
clear y_linha1
clear y_1
%% Estrat�gia 2 - Propulsor 35� + Propulsor 0� (n�o garante pouso suave em x)

%Altitude para ligar o propulsor
y_1 = ((y_linha_inicial^2) - (2*((E*cos(theta)/m) + g)*y_inicial))/(2*(E/m)*(1 - cos(theta)));
fprintf('Estrat�gia 2 - Propulsor 35� + Propulsor 0�: \n')
fprintf('A altitude de mudan�a de �ngulo � %.2f m.\n\n',y_1)
passo = floor((y_inicial-y_1)/0.001);
passo_1 = floor(y_1/0.001);

for i = 1:(passo+1)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha_inicial)^2) + (2*((E*cos(theta)/m) + g)*(y(i) - y_inicial)));
    x_linha1(i) = x_linha_inicial + ((y_linha1(i) - y_linha_inicial)*(E*sin(theta)/m)/((E*cos(theta)/m) + g));
end

for i = (passo+1):(passo + passo_1 + 2)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha1(passo + 1))^2) + (2*((E/m) + g)*(y(i) - y(passo + 1))));
    x_linha1(i) = x_linha1(passo+1) + ((y_linha1(i) - y_linha1(passo+1))*(E*sin(0)/m)/((E*cos(0)/m) + g));
end

figure(3)
plot(y_linha1,y)
grid
title('Estrat�gia 2 - Propulsor 35� + Propulsor 0�')
ylabel('Altitude (m)')
xlabel('Velocidade em y (m/s)')


figure(4)
plot(x_linha1,y)
grid
title('Estrat�gia 2 - Propulsor 35� + Propulsor 0�')
ylabel('Altitude (m)')
xlabel('Velocidade em x (m/s)')


%% Estrat�gia 3 - Propulsor a 35� + Propulsor a -35� + Propulsor a 0� (garante pouso suave em x e y)
y_2 = y_1;
v_1 = (y_linha_inicial + y_linha1(passo + 1))/2;
y_1 = (-(y_linha_inicial^2) + ((v_1)^2) + (2*((E*cos(theta)/m) + g)*y_inicial))/(2*((E*cos(theta)/m) + g));

clear y
clear y_linha1

fprintf('Estrat�gia 3 - Propulsor 35� + Propulsor -35� + Propulsor 0�: \n')
fprintf('A altitude da 1� mudan�a de �ngulo � %.2f m.\n',y_1)
fprintf('A altitude da 2� mudan�a de �ngulo � %.2f m.\n\n',y_2)
passo1 = floor((y_inicial - y_1)/0.001);
passo2 = floor((y_1 - y_2)/0.001);
passo3 = floor(y_2/0.001);

for i = 1:(passo1+1)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha_inicial)^2) + (2*((E*cos(theta)/m) + g)*(y(i) - y_inicial)));
    x_linha1(i) = x_linha_inicial + (((y_linha1(i) - y_linha_inicial)*(E*sin(theta)/m))/((E*cos(theta)/m) + g));
end

for i = (passo1 + 1):(passo1 + passo2 + 2)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha1(passo1 + 1))^2) + (2*((E*cos(-theta)/m) + g)*(y(i) - y(passo1 + 1))));
    x_linha1(i) = x_linha1(passo1 + 1) + ((y_linha1(i) - y_linha1(passo1 + 1))*(E*sin(-theta)/m)/((E*cos(-theta)/m) + g));
end

for i = (passo1 + passo2 + 2):(passo1 + passo2 + passo3 + 3)
    y(i) = y_inicial - (0.001*(i-1));
    y_linha1(i) = -sqrt(((y_linha1(passo1 + passo2 + 2))^2) + (2*((E/m) + g)*(y(i) - y(passo1 + passo2 + 2))));
    x_linha1(i) = x_linha1(passo1 + passo2 + 2) + ((y_linha1(i) - y_linha1(passo1 + passo2 + 2))*(E*sin(0)/m)/((E*cos(0)/m) + g));
end

figure(5)
plot(y_linha1,y)
grid
title('Estrat�gia 3 - Propulsor a 35� + Propulsor a -35� + Propulsor a 0�')
ylabel('Altitude (m)')
xlabel('Velocidade em y (m/s)')

figure(6)
plot(x_linha1,y)
grid
title('Estrat�gia 3 - Propulsor a 35� + Propulsor a -35� + Propulsor a 0�')
ylabel('Altitude (m)')
xlabel('Velocidade em x (m/s)')


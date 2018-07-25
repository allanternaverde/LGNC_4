% alunagem
clear all;
close all;
clc;
%% dados de entrada
% tracao do jato
E = 2.5;
% gravidade lunar
g = -1.5;
% massa do veiculo espacial
M = 1;
% condicoes de velocidade inicial para o grafico
v0 = [-20:0.1:0];
% valores de theta(inclinacao da propulsao) possíveis
theta = deg2rad([0 35]);
% vetore com as aceleracoes permitidas dado os angulos theta
a = E*cos(theta)/M+g;

% condicoes de pouso suave
v_f = 0;
y_f = 0;

% condicoes inicias de velocidade e altura
v_0 = -14;
y_0 = 150;


%% curvas de pouso suave
% torricelli com condiçoes finais zeros e jato a 0 graus
y0(:,1) = v0.^2/(2*a(1));
% torricellli com condições finais zeros e jato a 35 graus
y0(:,2) = v0.^2/(2*a(2));

% figura com condicoes de pouso suave
figure();
hold on;

% plot das curvas de pouso suave
plot(v0,y0,'Linewidth', 2)
leg = legend({'E à 0°','E à \pm35°'});
leg = leg.String;

%% condicoes iniciais do problema
%plot das condicoes inicias
scatter(v_0,y_0,'filled','Linewidth', 2, 'MarkerEdgeColor','k');
leg = [leg,'cond. iniciais (v_0=-14m/s; y_0=150m)'];
legend(leg);

%% condicao de chavemento queda-livre -> 0°
% aceleracao antes do chaveamento
a1 = g;
% aceleracao depois do chaveamento
a2 =  E*cos(deg2rad(0))/M+g;

% calcula o ponto de chaveamento
vals = [v_0, y_0, v_f, y_f, a1, a2];
[y_1, v_1, ltx] = ptoChaveamento(vals);

% plota o pronto de chaveamento
scatter(v_1(1),y_1(1),'filled','Linewidth', 2, 'MarkerEdgeColor','k');
leg = [leg,'queda-livre -> 0°'];

%% condicao de chavemento 0° -> 35°
% aceleracao antes do chaveamento
a1 =  E*cos(deg2rad(0))/M+g;
% aceleracao depois do chaveamento
a2 = E*cos(deg2rad(35))/M+g;

% calcula o ponto de chaveamento
vals = [v_0, y_0, v_f, y_f, a1, a2];
[y_2, v_2, ltx] = ptoChaveamento(vals);

% plota ponto de chaveamento
scatter(v_2(1),y_2(1),'filled','Linewidth', 2, 'MarkerEdgeColor','k');
leg = [leg,'0° -> 35°'];

%% condicao de chavemento 35° -> 0°
% aceleracao antes do chaveamento
a1 = E*cos(deg2rad(35))/M+g;
% aceleracao depois do chaveamento
a2 =  E*cos(deg2rad(0))/M+g;

% calcula o ponto de chaveamento
vals = [v_0, y_0, v_f, y_f, a1, a2];

% plota ponto de chaveamento
[y_3, v_3, ltx] = ptoChaveamento(vals);
scatter(v_3(1),y_3(1),'filled','Linewidth', 2, 'MarkerEdgeColor','k');
leg = [leg,'35° -> 0°'];

%% area de pouso possível
hf = fill_between(v0,y0(:,1),y0(:,2));
hf.FaceColor = [203 255 164]/255;
zero = zeros(length(y0(:,1)),1);

%% area de colisao
hf = fill_between(v0,zero,y0(:,2));
hf.FaceColor = [255 164 164]/255;

%% area de nao retorno porem pouco possível
ums = ones(length(y0(:,1)),1);
hf = fill_between(v0,y0(:,1),400*ums);
hf.FaceColor = [194 200 255]/255;

leg = ['pouso com queda livre inicial', 'colisão','pouso com acionamento inicial',leg];
legend(leg);
xlabel('v_0');
ylabel('y_0');

%% tebela com pontos de chaveamento
T = table({'queda livre';'colisão';'colisão'}, ...
    { ['(',num2str(v_1(1)),'m/s; ',num2str(y_1(1)),'m)']; 'não pousa' ; ['(',num2str(v_3(1)),'m/s; ',num2str(y_3(1)),'m)'] },...
    { 'colisão'; ['(',num2str(v_2(1)),'m/s; ',num2str(y_2(1)),'m)']; 'colisão' });
T.Properties.VariableNames = {'zero' 'Ecos0' 'Ecos35'};
T.Properties.RowNames = {'zero' 'Ecos0' 'Ecos35'};
writetable(T);
disp('Pontos de chaveamento:');
disp(T);


%% calculo do tempo até o primeiro chaveamento
a1 =  E*cos(deg2rad(0))/M+g;
a2 = E*cos(deg2rad(35))/M+g;

t(1) = (v_1(1) -v_0)/g;
t(2) = (v_2(1) -v_0)/a1;
t(3) = (v_3(1) -v_0)/a2;

%% tabelo com tempos ate primeiro chaveamento
T_t = table({'queda livre';'colisão';'colisão'}, ...
    { [num2str(t(1)),'s']; 'não pousa' ; [num2str(t(3)),'s'] },...
    { 'colisão'; [num2str(t(2)),'s']; 'colisão' });
T_t.Properties.VariableNames = {'zero'  'Ecos0' 'Ecos35'};
T_t.Properties.RowNames = {'zero'  'Ecos0' 'Ecos35'};
writetable(T_t);
disp('Tempos até o primeiro chaveamento:');
disp(T_t);

%% calculo do tempo desde o chaveamento até o pouso
t_2(1) = (0 - v_1(1))/a1;
t_2(2) = (0 - v_2(1))/a2;
t_2(3) = (0 - v_3(1))/a1;

T_t_2 = table({'queda livre';'colisão';'colisão'}, ...
    { [num2str(t_2(1)),'s']; 'não pousa' ; [num2str(t_2(3)),'s'] },...
    { 'colisão'; [num2str(t_2(2)),'s']; 'colisão' });
T_t_2.Properties.VariableNames = {'zero'  'Ecos0' 'Ecos35'};
T_t_2.Properties.RowNames = {'zero'  'Ecos0' 'Ecos35'};
writetable(T_t_2);


%% calculo do tempo total até o pouso
t_t(1) = t(1) + t_2(1);
t_t(2) = t(2) + t_2(2);
t_t(3) = t(3) + t_2(3);

T_t_t = table({'queda livre';'colisão';'colisão'}, ...
    { [num2str(t_t(1)),'s']; 'não pousa' ; [num2str(t_t(3)),'s'] },...
    { 'colisão'; [num2str(t_t(2)),'s']; 'colisão' });
T_t_t.Properties.VariableNames = {'zero'  'Ecos0' 'Ecos35'};
T_t_t.Properties.RowNames = {'zero'  'Ecos0' 'Ecos35'};
writetable(T_t_t);
disp('Tempo total até o pouso:');
disp(T_t_t);

%% plot das trajetórias
figure();
hold on;
grid minor;
leg = {};

%% pontos de chaveamento
scatter(v_0,y_0,'filled','Linewidth', 2, 'MarkerEdgeColor','k');
scatter(v_1(1),y_1(1),'filled','Linewidth', 2, 'MarkerEdgeColor','k');
scatter(v_2(1),y_2(1),'filled','Linewidth', 2, 'MarkerEdgeColor','k');
scatter(v_3(1),y_3(1),'filled','Linewidth', 2, 'MarkerEdgeColor','k');
plot(v0,y0,'Linewidth', 0.1)

%% trajetoria em queda livre a partir da posicao inicia
traj = trajetoria(g, v_0, y_0, t(1));
plot(traj.v, traj.y, 'Linewidth', 2, 'LineStyle', '-.', 'Color','g');

%% trajetoria com acionamento a 0 graus apos queda livre
traj = trajetoria(a1, v_1(1), y_1(1), t_2(1));
ln1 = plot(traj.v, traj.y, 'Linewidth', 2, 'LineStyle', '-.', 'Color','g');

%% trajetoria com acionamento a 0 graus a partir da posição inicial
traj = trajetoria(a1, v_0, y_0, t(2));
plot(traj.v, traj.y, 'Linewidth', 2, 'LineStyle', '-.', 'Color','m');

%% trajetoria com acionamento a 35 graus apos jato a 0 graus
traj = trajetoria(a2, v_2(1), y_2(1), t_2(2));
ln2 = plot(traj.v, traj.y, 'Linewidth', 2, 'LineStyle', '-.', 'Color','m');

%% trajetoria com acionamento a 35 graus a partir da posição inicial
traj = trajetoria(a2, v_0, y_0, t(3));
plot(traj.v, traj.y, 'Linewidth', 2, 'LineStyle', '-.', 'Color','b');

%% trajetoria com acionamento a 0 graus apos jato a 35 graus
traj = trajetoria(a1, v_3(1), y_3(1), t_2(3));
ln3 = plot(traj.v, traj.y, 'Linewidth', 2, 'LineStyle', '-.', 'Color','b');


lines = [ln1 ln2 ln3];
legds = {'1.3832s + 16.0748s = 17.458s',...
        '2.7738s + 20.4902s = 23.264s',...
        '10.2852s + 8.365s = 18.6501s'};
    
legend(lines, legds);
xlabel('v_0');
ylabel('y_0');


%% Estrat�gia 3 - Propulsor a 35� + Propulsor a -35� + Propulsor a 0� (garante pouso suave em x e y)
m = M;
theta = 35*pi/180;
y_inicial = 150;
y_linha_inicial = -14;
x_inicial = 0;
x_linha_inicial = 0;

y_1 = ((y_linha_inicial^2) - (2*((E*cos(theta)/m) + g)*y_inicial))/(2*(E/m)*(1 - cos(theta)));
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

figure();
hold on;
plot(x_linha1,y,'Linewidth', 2, 'LineStyle', '-.', 'Color','b')
grid minor;
ylabel('y [m]')
xlabel('v_x[m/s]')

y_2 = y_1;
v_1 = (y_linha_inicial + y_linha1(passo + 1))/2;
y_1 = (-(y_linha_inicial^2) + ((v_1)^2) + (2*((E*cos(theta)/m) + g)*y_inicial))/(2*((E*cos(theta)/m) + g));

clear y
clear y_linha1

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

figure();
hold on;
plot(x_linha1,y,'Linewidth', 2, 'LineStyle', '-.', 'Color','b')
grid minor;
ylabel('y [m]')
xlabel('v_x [m/s]')


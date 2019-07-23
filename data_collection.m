set(0,'DefaultFigureVisible', 'on')
dfdir = 'C:\Users\kxxs\Desktop\Energy-Hub\Load\';
rand('seed',1);

Time = 1:24;
            %1    2      3    4      5    6     7     8
Price_E = [0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.64, 1.12,...
           1.12, 1.12, 1.12, 1.12, 0.64, 0.64, 0.64, 0.64,...
           1.12, 1.12, 1.12, 1.12, 0.64, 0.64, 0.64, 0.31];

          %1    2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
          %17   18  19  20  21  22  23  24
Solar  = [ 0,   0,  0,  0,  0,  3,  5,  16, 26, 34, 38, 40, 42, 40, 36, 26,...
           12,  4,  1,  0,  0,  0,  0,  0;
          ]*1;

%% Summer Load
Temp_Outside_Summer  = [ 29.5, 29.1, 29, 28.4, 28.2, 28.3, 28.3, 29, 29.9, 30.5, 31.6, 32.3, 33.1, 33.6, 33.8, 33.9,...
           34, 33.2, 32.6, 31.9, 31.3, 30.3, 29.6, 28.3;]; % 2005.7.16 Nanjing temperature
%        40,   39, 39.6, 40, 40, 39.2, 40.2, 42,48, 60, 76, 100, 120, 108, 112, 108,110, 92, 80,...
%           88,   90, 84,   62, 48;  tongli
Demand_Summer = [45,   47, 46, 44, 45, 53, 56, 85, 104, 110, 120, 115, 92, 87, 89, 85,...
          86,   89, 92, 91, 79, 70, 65, 55;    % W (electricity)
          17,   16, 16, 16, 16, 18, 20, 25, 25, 27, 30, 37, 35, 34, 33, 32....
          30,   29, 29, 29, 27, 23, 23, 24;   % R (cold)
          30,   30, 30, 30, 32, 32, 36, 40, 42, 35, 32, 28, 24, 20, 21, 23,...
          28,   32, 25, 26, 27, 30, 32, 35;   % Q (head)
          ]; 
Demand_Summer(2,:) = (Temp_Outside_Summer - 25) * 8;  % cooling demand estimation
Demand_Summer(3,:) = max(30 - Temp_Outside_Summer,0) * 4 + 15; % heat demand estimation

%% Winter Load
Temp_Outside_Winter  = [1.2,1.2,1.5,1.4,1.5,1.5,1.3,1.1,1.1,1.3,1.7,2.3...
                        2.8,2.9,3.3,3.1,2.8,2.3,2.2,2.0,1.8,1.4,1.2,1.2]; % 2005.1.15 Nanjing temperature
%        40,   39, 39.6, 40, 40, 39.2, 40.2, 42,48, 60, 76, 100, 120, 108, 112, 108,110, 92, 80,...
%           88,   90, 84,   62, 48;  tongli
Demand_Winter = [45,   47, 46, 44, 45, 53, 56, 85, 104, 110, 108, 115, 92, 87, 89, 85,...
          86,   89, 92, 91, 79, 70, 65, 55;   % W (electricity)
          17,   16, 16, 16, 16, 18, 20, 25, 25, 27, 30, 37, 35, 34, 33, 32....
          30,   29, 29, 29, 27, 23, 23, 24;   % R (cold)
          30,   30, 30, 30, 32, 32, 36, 40, 42, 35, 32, 28, 24, 20, 21, 23,...
          28,   32, 25, 26, 27, 30, 32, 35;   % Q (head)
          ]; 
Demand_Winter(1,:) = (Demand_Winter(1,:) - mean(Demand_Winter(1,:)))/1.5 + mean(Demand_Winter(1,:)); 
Demand_Winter(2,:) = (rand(1,24) - 0.5) * 5 + 5;  % cooling demand estimation
Demand_Winter(3,:) =  (6 - Temp_Outside_Winter)* 8 + 15; % heat demand estimation

%%
Temp_Outside_Spring = [3.5,3,2.1,1.5,1.1,0.6,0.7,2.4,5.8,8,10.5,12.2,13.7,14.8,...
                15.2,15.7,15.8,15.9,14.1,12.9,12,11.5,11.3,11];

Demand_Spring  = [45,   47, 46, 44, 45, 53, 56, 85, 104, 110, 108, 115, 92, 87, 89, 85,...
          86,   89, 92, 91, 79, 70, 65, 55;   % W (electricity)
          17,   16, 16, 16, 16, 18, 20, 25, 25, 27, 30, 37, 35, 34, 33, 32....
          30,   29, 29, 29, 27, 23, 23, 24;   % R (cold)
          30,   30, 30, 30, 32, 32, 36, 40, 42, 35, 32, 28, 24, 20, 21, 23,...
          28,   32, 25, 26, 27, 30, 32, 35;   % Q (head)
          ];  
Demand_Spring(1,:) = Demand_Winter(1,:) * 2/3 + Demand_Summer(1,:) *1/3;
Demand_Spring(2,:) = max(Temp_Outside_Spring - 12,0) / 2 + 5;  % cooling demand estimation
Demand_Spring(3,:) = max(14 - Temp_Outside_Spring,0) + 12;  % heat demand estimation
%% plot
figure
% subplot(1,3,1)
% plot(Time,Solar,'-*');
% hold on
plot(Time,Demand_Summer(1,:),'-.');
hold on
plot(Time,Demand_Summer(2,:));
hold on
plot(Time,Demand_Summer(3,:));
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
% 'solar',
legend('elec','cold','heat')
title('Original Load (Summer)')
saveas(gcf, [dfdir,'load_summer','.jpg'])


% subplot(1,3,2)
figure
% plot(Time,Solar,'-*');
% hold on
plot(Time,Demand_Winter(1,:),'-.');
hold on
plot(Time,Demand_Winter(2,:));
hold on
plot(Time,Demand_Winter(3,:));
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
% 'solar',
legend('elec','cold','heat')
title('Original Load (Winter)')
saveas(gcf, [dfdir,'load_winter','.jpg'])

% subplot(1,3,3)
figure
% plot(Time,Solar,'-*');
% hold on
plot(Time,Demand_Spring(1,:),'-.');
hold on
plot(Time,Demand_Spring(2,:));
hold on
plot(Time,Demand_Spring(3,:));
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
% 'solar',
legend('elec','cold','heat')
title('Original Load (Spring)')
saveas(gcf, [dfdir,'load_spring','.jpg'])

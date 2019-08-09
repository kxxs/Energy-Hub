% version 5 : potential assessment
% Demand of W,R,Q, 按照BT由小到大排

% renewables = importdata('renewable.dat');
% Wind_Summer = renewables.Wind_summer_unit;
% Wind_Winter = renewables.Wind_winter_unit;

fdir = 'C:\Users\kxxs\Desktop\Energy-Hub\Summer\';
loaddir = 'C:\Users\kxxs\Desktop\Energy-Hub\Load\';

% data_collection;
Demand = Demand_Summer;
% load([loaddir,'Demand_Spring.mat']);

mysolver = 'cplex';

Gas_Coef = 1;
E = zeros(1,25);  % storage amount
Storage_Cap = [40,40,40];  % storage capacity
Storage_Pm = [20,20,20];   % max power of storage

Gas_In = zeros(1,24);
Elec_In = zeros(1,24);
%%
%%%%%%%%%%%%%%%% Branch Matrix (based on the paper) %%%%%%%%%%%%%%%%%%%%%%%%%%
% BT(Branch Type): 1-W;  2-R(cold);  3-Q;  4-F(gas) 5-solar 6-wind
% s: start node;
% t: end node;
No = 1; BT = 2; s = 3; t = 4; cap = 5;
Branch = [
    % BT(Branch Type): 1-W;  2-R(cold);  3-Q;  4-F(gas) 5-solar 6-wind
    % No.  BT   s   t   capacity
       1    1  -1   0   120;
       2    1   -1  3   10;
       3    1   -1  4   10;
       4    4   -1  1   55 * Gas_Coef;
       5    4   -1  2   20 * Gas_Coef;
       6    1   1   3   10;
       7    1   1   4   0;
       8    1   1   0   80;
       9    1   1   7   80;
       10   3   1   0   80;
       11   3   1   8   40;
       12   3   2   6   0;
       13   3   2   0   90;
       14   2   3   0   90;
       15   2   3   5   90;
       16   2   5   0   90;
       17   3   4   0   90;
       18   3   4   6   90;
       19   3   6   0   100;
       20   1   7   0   100;
       21   2   8   0   100;
       22   1   -1  7   20;
       23   3   1   6   20;
       24   2   5   inf 90;    % this branch is for △E  (dE), such branches should be at last
       25   3   6   inf 90;
       26   1   7   inf 90;
       ];
Branch_Num = max(Branch(:,1));
Input_Num = length(Branch(Branch(:,s)==-1,1)); % input num of the energy hub
Output_Num = length(Branch(Branch(:,t)==0,1));

Input_BT = unique(Branch(Branch(:,s)==-1,BT));
Output_BT = unique(Branch(Branch(:,t)==0,BT));

V = sdpvar(length(Branch(Branch(:,t)~=inf,1)),24);
dE = sdpvar(length(Branch(Branch(:,t)==inf,1)),24);
V_In = sdpvar(length(Input_BT),24);
V_Out = sdpvar(length(Output_BT),24); 
Cutdown = sdpvar(3,24);
Shift = sdpvar(3,24);
SolarUsed = sdpvar(1,24);
%%
%%%%%%%%%%%%%% Node Matrix %%%%%%%%%%%%%%%%%
NT = 2; p1 = 3; p2 = 4;
Node = [
% NT = 1: input coef = eta, output coef = 1 in Z matrix
% NT = 2: input coef = eta1, output coef = eta2, storage coef = 1 (storage node)
%  No.  NT   p1      p2
   -1   1    0       0;
   0    1    0       0;
   1    1    0.35    0.45;  %% 哪个branch的index小，该branch的能量转换效率就写在前面
   2    1    0.85    0;
   3    1    2.90    0;
   4    1    2.30    0;
   5    2    0.90    0.95;
   6    2    0.90    0.95;
   7    2    0.90    0.95;
   8    1    0.90    0;
       ];  

Cons = [sum(Shift(1,:)) == 0;   % 削峰填谷约束
		sum(Shift(2,:)) == 0;
		sum(Shift(3,:)) == 0;
        sum(dE(1,:)) == 0;  %  储能平衡约束
        sum(dE(2,:)) == 0;
        sum(dE(3,:)) == 0;];
   %%
%%%%%%%%%%%%%%%%%% Main %%%%%%%%%%%%%%%%%%%%%%%%%
for hour = 1:24
    Node_Num = max(Node(:,1));
    Z = [];
    Zs = []; % storage matrix
    for i = 1:Node_Num
        if(Node(Node(:,1)==i,NT)==2)
            Is_Storage = 1;
            Storage_Dim = 1;
        else
            Is_Storage = 0;
            Storage_Dim = 0;
        end
        % find branches pointing into and from current node 
        In_B = Branch(Branch(:,t)==i,:); In_B_backup = In_B;
        Out_B = Branch(Branch(:,s)==i,:); 
        if(Is_Storage == 1)
            for Out_Idx = 1:length(Out_B(:,1))
                if(Out_B(Out_Idx,t)==inf)
                    Storage_B = Out_B(Out_Idx,:);
                    Out_B(Out_Idx,:) = [];
                end
            end
        end
        Out_B_backup = Out_B;
        % find num of different energy type
        In_Dim = length(unique(In_B(:,BT)));
        Out_Dim = length(unique(Out_B(:,BT)));
        % construct H and A matrix
        if(Is_Storage == 1)
            H = zeros(1, In_Dim + Out_Dim);
        else
            H = zeros(In_Dim * Out_Dim, In_Dim + Out_Dim + Storage_Dim);
        end
        A = zeros(In_Dim + Out_Dim + Storage_Dim , Branch_Num);
        % fill the elements of A
        for A_idx = 1:length(A(:,No))
            if(isempty(In_B)==0) % fill input energy into matrix A first
                % ET: energy type
                ET = In_B(1,BT);
                A(A_idx, In_B(In_B(:,BT)==ET,No)) = 1;
                In_B(In_B(:,BT)==ET,:) = [];
            elseif(isempty(Out_B)==0)  % then fill output energy into matrix A
                ET = Out_B(1,BT);
                A(A_idx, Out_B(Out_B(:,BT)==ET,No)) = -1;
                Out_B(Out_B(:,BT)==ET,:) = [];
            else % fill storage energy into A
                A(A_idx, Storage_B(1,No)) = -1;
            end
        end
        % fill the elements of H (preparation)
        In_B = In_B_backup;
        Out_B = Out_B_backup;
        % fill the elements of H 
        if(Is_Storage == 1) % storage node
            H(1,1) = Node(Node(:,1)==i,p1);
            H(1,2) = 1/Node(Node(:,1)==i,p2);
            H(1,3) = 1;
        else % not storage node
            % list of in/out energy type
            In_temp=In_B(:,BT); In_temp = In_temp(end:-1:1);
            Out_temp=Out_B(:,BT); Out_temp = Out_temp(end:-1:1);
            [In_List, In_rank] = unique(In_temp); 
            [Out_List, Out_rank] = unique(Out_temp);  
            In_List = In_temp(sort(In_rank)); In_List = In_List(end:-1:1); % keep the index order, Branch index小的能量形式在前，index 大的在后
            Out_List = Out_temp(sort(Out_rank)); Out_List = Out_List(end:-1:1);
            H_Idx = 1; Param_Idx = p1;
            for j = 1:length(In_List)
                for k = 1:length(Out_List)
                    H(H_Idx,j) = Node(Node(:,1)==i,Param_Idx);
                    H(H_Idx,k + length(In_List)) = 1;
                    H_Idx = H_Idx + 1;
                    Param_Idx = Param_Idx + 1;
                    %%% 这里考虑H矩阵每一行都是一种输入(系数为eta)、一种输出（系数为1），其余都是0，不知道有没有问题
                    %%% 然后把j种输入,k种输出的j*k个组合枚举了一遍。
                end
            end
        end
        Z = [Z;H*A];
    end
    %% 
    % Construct X and Y

    X = zeros(length(Input_BT), Branch_Num);
    Y = zeros(length(Output_BT), Branch_Num);

    for i = 1:length(Input_BT)
        X(i,intersect(Branch(Branch(:,s)==-1),Branch(Branch(:,BT)== Input_BT(i)))) = 1;
    end

    for i = 1:length(Output_BT)
        Y(i,intersect(Branch(Branch(:,t)==0),Branch(Branch(:,BT)== Output_BT(i)))) = 1;
    end

    %%%%%%%%%%%%%%%%%%% solve with yalmip %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Cons = [ Cons;
            [X;Y;Z]*[V(:,hour);dE(:,hour)] == [V_In(:,hour);V_Out(:,hour);zeros(length(Z(:,1)),1)];
            SolarUsed(1,hour) <= Solar(hour);   % 光伏出力约束
            SolarUsed(1,hour) >= 0.5*Solar(hour); % 要求光伏利用率>=0.6
            SolarUsed(1,hour) <= V_In(1,hour);
            dE(:,hour) <= [Storage_Pm(1);Storage_Pm(2);Storage_Pm(3)];   
            dE(:,hour) >= [-Storage_Pm(1);-Storage_Pm(2);-Storage_Pm(3)];  % 充放电功率约束
            sum(dE(1,1:hour)) <= Storage_Cap(1);  
            sum(dE(2,1:hour)) <= Storage_Cap(2); 
            sum(dE(3,1:hour)) <= Storage_Cap(3);  % 储能容量约束
            sum(dE(1,1:hour)) >= -Storage_Cap(1)*0.8;
            sum(dE(2,1:hour)) >= -Storage_Cap(2)*0.8;  
            sum(dE(3,1:hour)) >= -Storage_Cap(3)*0.8;
 			Cutdown(:,hour) <= 0.1*Demand(:,hour);  % 负荷削减约束
            Cutdown >=0;
            Shift(:,hour) <= 10;
 			Shift(:,hour) <= 0.3*Demand(:,hour);    % 负荷转移约束
 			Shift(:,hour) >= -0.1*Demand(:,hour);
            V_Out(:,hour) == (Demand(:,hour) - Cutdown(:,hour) + Shift(:,hour));
            V_In(:,hour) >=0; V_Out(:,hour) >=0; V(:,hour)>=0;	
            [V(:,hour);dE(:,hour)] <= Branch(:,cap);
             V_In(1,:) - SolarUsed(1,:) <= 80; % 供应量约束
             V_In(2,:) <= 50 * Gas_Coef;
            ];
        
         if(hour > 1)
         Cons = [Cons;  % Change rate constraint
              -20 <= V(:,hour) -   V(:,hour-1) <= 20;
              -15 <= V(8,hour) -   V(8,hour-1) <= 15;   % AB
              -10 <= V(2,hour) -   V(2,hour-1) + V(6,hour) -   V(6,hour-1) <= 10;  % EC
              -10 <= V(3,hour) -   V(3,hour-1) + V(7,hour) -   V(7,hour-1) <= 10;  % EH
             ];
         end
% final2.CHP = double(V(4,:));
% final2.GB = double(V(5,:));
% final2.EC = double(V(2,:) + V(6,:));
% final2.EH = double(V(3,:) + V(7,:));
% final2.AB = double(V(8,:));
end
z = binvar(3,24);
Shift_Abs = sdpvar(3,24);
U = 1000;


Cons1 = [Cons;
        0 <= Shift_Abs - Shift <= 2*U*z; U*(1-z) >= Shift;
        0 <= Shift_Abs + Shift <= 2*U*(1-z); -U*z <= Shift;
        -U <= Shift <= U;];

Cost = (V_In(1,:) - SolarUsed(1,:))*Price_E' + sum(2.85/10*V_In(2,:));
ops = sdpsettings('solver',mysolver,'verbose',1);
solvesdp(Cons1,Cost,ops);

final.Demand = Demand;
final.ElecInput = double(V_In(1,:));
final.GasInput = double(V_In(2,:));
final.ElecOutput = double(V_Out(1,:));
final.ColdOutput = double(V_Out(2,:));
final.HeatOutput = double(V_Out(3,:));
final.ElecShift = double(Shift(1,:));
final.ColdShift = double(Shift(2,:));
final.HeatShift = double(Shift(3,:));
final.ElecShift_Abs = double(Shift_Abs(1,:));
final.ColdShift_Abs = double(Shift_Abs(2,:));
final.HeatShift_Abs = double(Shift_Abs(3,:));
final.ElecCutdown = double(Cutdown(1,:));
final.ColdCutdown = double(Cutdown(2,:));
final.HeatCutdown = double(Cutdown(3,:));
final.Cost = double((V_In(1,:) - SolarUsed(1,:))*Price_E' + sum(2.05/10*V_In(2,:))/2);
final.TotalCost = double(Cost);
final.Storage = double(dE);
final.SolarUsed = double(SolarUsed);
final.CHP = double(V(4,:));
final.GB = double(V(5,:));
final.EC = double(V(2,:) + V(6,:));
final.EH = double(V(3,:) + V(7,:));
final.AB = double(V(8,:));
final.V = double(V);

% Total Reduction of Electricity Demand during peak periods after response 
final.OutputElecCut = [Demand(1,8:12),Demand(1,17:19)] - [double(V_Out(1,8:12)),double(V_Out(1,17:19))];

%%
ops = sdpsettings('solver',mysolver,'verbose',1);
Cons2 = [Cons; Shift == 0; Cutdown ==0;];
solvesdp(Cons2,Cost,ops);

final2.Demand = Demand;
final2.ElecInput = double(V_In(1,:));
final2.GasInput = double(V_In(2,:));
final2.ElecOutput = double(V_Out(1,:));
final2.ColdOutput = double(V_Out(2,:));
final2.HeatOutput = double(V_Out(3,:));
final2.ElecShift = double(Shift(1,:));
final2.ColdShift = double(Shift(2,:));
final2.HeatShift = double(Shift(3,:));
final2.ElecCutdown = double(Cutdown(1,:));
final2.ColdCutdown = double(Cutdown(2,:));
final2.HeatCutdown = double(Cutdown(3,:));
final2.Cost = double(Cost);
final2.TotalCost = double(Cost2);
final2.Storage = double(dE);
final2.SolarUsed = double(SolarUsed);
final2.CHP = double(V(4,:));
final2.GB = double(V(5,:));
final2.EC = double(V(2,:) + V(6,:));
final2.EH = double(V(3,:) + V(7,:));
final2.AB = double(V(8,:));

% Total Reduction of Supplied Electricity during peak periods after response 
final.InputElecCut = [final2.ElecInput(8:12),final2.ElecInput(17:19)]...
                          - [final.ElecInput(8:12),final.ElecInput(17:19)];

%% plot the optimal input for each hour
set(0,'DefaultFigureVisible', 'on')

figure
plot(Time,final.ElecInput); hold on;
plot(Time,final.GasInput); hold on;
plot(Time,final.SolarUsed); hold on;
plot(Time,Solar);
axis([0 25 0 140])
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('Elec Input','Gas Input','Solar Used','Solar Total')
title('Input graph')
saveas(gcf, [fdir,'input_withIDR','.jpg'])

figure
plot(Time,final2.ElecInput); hold on;
plot(Time,final2.GasInput); hold on;
plot(Time,final2.SolarUsed); hold on;
plot(Time,Solar);
axis([0 25 0 140])
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('Elec Input','Gas Input','Solar Used','Solar Total')
title('Input graph without IDR')
saveas(gcf, [fdir,'input_noIDR','.jpg'])

figure
plot(Time, Demand(1,:)); hold on; plot(Time, final.ElecOutput); hold on;
plot(Time, Demand(2,:)); hold on; plot(Time, final.ColdOutput); hold on;
plot(Time, Demand(3,:)); hold on; plot(Time, final.HeatOutput);
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('ElecOrig','ElecResp','ColdOrig','ColdResp','HeatOrig','HeatResp','Location','NorthWest')
title('Output Graph')
saveas(gcf, [fdir,'output','.jpg'])

% figure
% figure
% plot(Time, final2.ElecOutput); hold on; plot(Time, final.ElecOutput); hold on;
% plot(Time, final2.ColdOutput); hold on; plot(Time, final.ColdOutput); hold on;
% plot(Time, final2.HeatOutput); hold on; plot(Time, final.HeatOutput);
% legend('ElecOrig','ElecResp','ColdOrig','ColdResp','HeatOrig','HeatResp','Location','NorthWest')
% title('Output Graph')


figure
plot(Time,final.ElecShift); hold on;
plot(Time,final.ColdShift); hold on; 
plot(Time,final.HeatShift);
axis([0 25 -30 30])
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('ElecShift','ColdShift','HeatShift')
title('Shift Amount')
saveas(gcf, [fdir,'shift','.jpg'])

figure
plot(Time,final.ElecCutdown); hold on;
plot(Time,final.ColdCutdown); hold on;
plot(Time,final.HeatCutdown);
axis([0 25 0 30])
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('ElecCutdown','ColdCutdown','HeatCutdown')
title('Cutdown Amount')
saveas(gcf, [fdir,'cutdown','.jpg'])

figure
plot(Time,final.CHP); hold on;
plot(Time,final.GB); hold on;
plot(Time,final.EC); hold on;
plot(Time,final.EH); hold on;
plot(Time,final.AB); hold on;
axis([0 25 0 70*Gas_Coef])
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('CHP','GB','EC','EH','AB')
title('Energy flow')
saveas(gcf, [fdir,'energy flow','.jpg'])

figure
plot(Time,final.Storage(1,:)); hold on;
plot(Time,final.Storage(2,:)); hold on;
plot(Time,final.Storage(3,:)); 
axis([0 25 -25 25])
plot([8 8], get(gca, 'YLim'), '--g')
plot([12 12], get(gca, 'YLim'), '--g')
plot([17 17], get(gca, 'YLim'), '--g')
plot([21 21], get(gca, 'YLim'), '--g')
legend('Electricity','Cold','Heat')
title('Energy Storage')
saveas(gcf, [fdir,'storage','.jpg'])

% Gain Coefficient = total reduction of elec supply / elec demand during peak periods
gain_coef = sum(final.InputElecCut)/sum(final.OutputElecCut)

% version 3 : demand response added
% Demand of W,R,Q, 按照BT由小到大排
E = zeros(1,25);  % storage amount
Time = 1:24;
          %1    2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
          %17   18  19  20  21  22  23  24
Solar  = [ 0,   0,  0,  0,  0,  3,  5,  16, 26, 34, 38, 40, 42, 40, 36, 26,...
           12,  4,  1,  0,  0,  0,  0,  0;
          ]*5;
Demand = [45,   47, 46, 44, 45, 53, 56, 77, 80, 92, 95, 94, 85, 83, 84, 82,...
          89,   101,110,115,110,90, 70, 50;   % W (electricity)
          17,   16, 16, 16, 16, 18, 20, 25, 25, 27, 30, 37, 35, 34, 33, 32....
          30,   29, 29, 29, 27, 23, 23, 24;   % R (cold)
          30,   30, 30, 30, 32, 32, 36, 40, 42, 35, 32, 28, 16, 17, 21, 23,...
          28,   32, 25, 26, 27, 30, 32, 35;   % Q (head)
          ]; 

figure(1)
plot(Time,Solar,'-*');
hold on
plot(Time,Demand(1,:),'-.');
hold on
plot(Time,Demand(2,:));
hold on
plot(Time,Demand(3,:));
legend('solar','elec','cold','heat')

Gas_In = zeros(1,24);
Elec_In = zeros(1,24);
%%
%%%%%%%%%%%%%%%% Branch Matrix (based on the paper) %%%%%%%%%%%%%%%%%%%%%%%%%%
% BT(Branch Type): 1-W;  2-R(cold);  3-Q;  4-F(gas) 5-solar 6-wind
% s: start node;
% t: end node;
No = 1; BT = 2; s = 3; t = 4; cap = 5;
% Branch = [
%     % No.  BT   s   t   capacity
%        1    4  -1   1   inf;
%        2    3   1   2   inf;
%        3    3   1   0   300;
%        4    1   1   0   inf;
%        5    2   2   0   inf;
%           ];
Branch = [
    % BT(Branch Type): 1-W;  2-R(cold);  3-Q;  4-F(gas) 5-solar 6-wind
    % No.  BT   s   t   capacity
       1    1  -1   0   inf;
       2    1  -1   2   inf;
       3    4  -1   1   inf;
       4    4  -1   3   inf;
       5    1   1   2   inf;
       6    1   1   0   inf;
       7    1   1   5   inf;
       8    3   1   0   inf;
       9    3   1   4   inf;
       10   3   3   4   inf;
       11   3   3   0   inf;
       12   2   2   0   inf;
       13   3   5   0   inf;
       14   2   4   0   inf;
       15   1   1   6   inf;
       16   3   6   0   inf;
       17   3   6   4   inf;
       18   3   6   inf inf;    % this branch is for △E  (dE), such branches should be at last
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
SolarUsed = sdpvar(1,24);
%%
%%%%%%%%%%%%%% Node Matrix %%%%%%%%%%%%%%%%%
NT = 2; p1 = 3; p2 = 4;
% Node = [
% %  No.  NT   p1      p2
%    -1  -1    0       0;
%    0    0    0       0;
%    1    3    0.25    0.50;  %% 哪个branch的index小，该branch的能量转换效率就写在前面
%    2    1    0.80    0;
%        ];
Node = [
% NT = 1: input coef = eta, output coef = 1 in Z matrix
%  No.  NT   p1      p2
   -1   1    0       0;
   0    1    0       0;
   1    1    0.50    0.25;  %% 哪个branch的index小，该branch的能量转换效率就写在前面
   2    1    1.20    0;
   3    1    0.80    0;
   4    1    0.90    0;
   5    1    0.30    0;
   6    2    0.90    0.90;
       ];  

Cons = [];
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
            SolarUsed(1,hour) >= 0.6*Solar(hour); % 要求光伏利用率>=0.6
            SolarUsed(1,hour) <= V_In(1,hour);
            dE(:,hour) <= 10;   dE(:,hour) >= -10;  % 充放电功率约束
            sum(dE(:,1:hour)) <= 50;    % 储能容量约束
            sum(dE(:,1:hour)) >= 0;
            V_Out(:,hour) >= Demand(:,hour);
            V_In(:,hour) >=0; V_Out(:,hour) >=0; V(:,hour)>=0;	
            [V(:,hour);dE(:,hour)] <= Branch(:,cap);
            ];
end
Cost = sum(1.2*(V_In(1,:) - SolarUsed(1,:)) + 2.05*V_In(2,:)); %Vin的次序也是按BT编号由小到大，如本例中1-W 4-Gas

ops = sdpsettings('solver','gurobi','verbose',1);
solvesdp(Cons,Cost,ops);

final.Demand = Demand;
final.ElecInput = double(V_In(1,:));
final.GasInput = double(V_In(2,:));
final.ElecOutput = double(V_Out(1,:));
final.ColdOutput = double(V_Out(2,:));
final.HeatOutput = double(V_Out(3,:));
final.Cost = double(Cost);
final.Storage = double(dE);
final.SolarUsed = double(SolarUsed);
for hour = 1:24
result(hour).V = double(V(:,hour));
end
%% plot the optimal input for each hour
figure
plot(Time,final.ElecInput); hold on;
plot(Time,final.GasInput); hold on;
plot(Time,final.SolarUsed); hold on;
plot(Time,Solar);
legend('Elec Input','Gas Input','Solar Used','Solar Total')
title('Input graph')

figure
plot(Time, Demand(1,:)); hold on; 
plot(Time, Demand(2,:)); hold on; 
plot(Time, Demand(3,:)); 
legend('ElecDemand','ColdDemand','HeatDemand')
title('Output Graph')



    

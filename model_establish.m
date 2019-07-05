Demand = [100;50;30]; % Demand of W,R,Q

%%%%%%%%%%%%%%%% Branch Matrix (based on the paper) %%%%%%%%%%%%%%%%%%%%%%%%%%
% BT(Branch Type): 1-W;  2-R(cold);  3-Q;  4-F(gas) 5-solar 6-wind
% s: start node;
% t: end node;
No = 1; BT = 2; s = 3; t = 4; cap = 5;
Branch = [
    % No.  BT   s   t   capacity
       1    4  -1   1   inf;
       2    3   1   2   inf;
       3    3   1   0   300;
       4    1   1   0   inf;
       5    2   2   0   inf;
          ];
Branch_Num = max(Branch(:,1));
Input_Num = length(Branch(Branch(:,s)==-1,1)); % input num of the energy hub
Output_Num = length(Branch(Branch(:,t)==0,1));

%%%%%%%%%%%%%% Node Matrix %%%%%%%%%%%%%%%%%
NT = 2; p1 = 3; p2 = 4;
Node = [
%  No.  NT   p1      p2
   -1  -1    0       0;
   0    0    0       0;
   1    3    0.25    0.50;  %% 哪个branch的index小，该branch的能量转换效率就写在前面
   2    1    0.80    0;
       ];
    
%%%%%%%%%%%%%%%%%% Main %%%%%%%%%%%%%%%%%%%%%%%%%
Node_Num = max(Node(:,1));
Z = [];
for i = 1:Node_Num
    % find branches pointing into and from current node
    In_B = Branch(Branch(:,t)==i,:); In_B_backup = In_B;
    Out_B = Branch(Branch(:,s)==i,:); Out_B_backup = Out_B;
    % find num of different energy type
    In_Dim = length(unique(In_B(:,BT)));
    Out_Dim = length(unique(Out_B(:,BT)));
    % construct H and A matrix
    H = zeros(In_Dim * Out_Dim, In_Dim + Out_Dim);
    A = zeros(In_Dim + Out_Dim, Branch_Num);
    % fill the elements of A
    for A_idx = 1:length(A(:,No))
        if(isempty(In_B)==0) % fill input energy into matrix A first
            % ET: energy type
            ET = In_B(1,BT);
            A(A_idx, In_B(In_B(:,BT)==ET,No)) = 1;
            In_B(In_B(:,BT)==ET,:) = [];
        else
            ET = Out_B(1,BT);
            A(A_idx, Out_B(Out_B(:,BT)==ET,No)) = -1;
            Out_B(Out_B(:,BT)==ET,:) = [];
        end
    end
    % fill the elements of H (preparation)
    In_B = In_B_backup;
    Out_B = Out_B_backup;
    % list of in/out energy type
    In_temp=In_B(:,BT); In_temp = In_temp(end:-1:1);
    Out_temp=Out_B(:,BT); Out_temp = Out_temp(end:-1:1);
    [In_List, In_rank] = unique(In_temp); 
    [Out_List, Out_rank] = unique(Out_temp);  
    In_List = In_List(sort(In_rank)); In_List = In_List(end:-1:1); % keep the index order, Branch index小的能量形式在前，index 大的在后
    Out_List = Out_List(sort(Out_rank)); Out_List = Out_List(end:-1:1);
%     In_Out_List = zeros(length(In_List) * length(Out_List), 2); % list all the in-out energy pairs for current node
    % fill the elements of H 
    H_Idx = 1; Param_Idx = p1;
    for j = 1:length(In_List)
        for k = 1:length(Out_List)
%             In_Out_List((j-1)*length(Out_List)+k,1) = In_List(j);
%             In_Out_List((j-1)*length(Out_List)+k,2) = Out_List(k);
%             Current_In =  In_List(j); % 当前的输入Hub的能量形式
%             Current_Out = Out_List(k); % 当前的从Hub输出的能量形式
            H(H_Idx,j) = Node(Node(:,1)==i,Param_Idx);
            H(H_Idx,k + length(In_List)) = 1;
            H_Idx = H_Idx + 1;
            Param_Idx = Param_Idx + 1;
            %%% 这里考虑H矩阵每一行都是一种输入(系数为eta)、一种输出（系数为1），其余都是0，不知道有没有问题
            %%% 然后把j种输入,k种输出的j*k个组合枚举了一遍。
        end
    end
    Z = [Z;H*A];
end

% Construct X and Y
X = zeros(Input_Num, Branch_Num);
Y = zeros(Output_Num, Branch_Num);

Input_BT = unique(Branch(Branch(:,s)==-1,BT));
Output_BT = unique(Branch(Branch(:,t)==0,BT));

for i = 1:Input_Num
    X(i,intersect(Branch(Branch(:,s)==-1),Branch(Branch(:,BT)== Input_BT(i)))) = 1;
end

for i = 1:Output_Num
    Y(i,intersect(Branch(Branch(:,t)==0),Branch(Branch(:,BT)== Output_BT(i)))) = 1;
end

V = sdpvar(length(Branch),1);
V_In = sdpvar(Input_Num,1);
V_Out = sdpvar(Output_Num,1);

%%%%%%%%%%%%%%%%%%% solve with yalmip %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Cost = sum(100*V_In);
Cons = [
        [X;Y;Z]*V == [V_In;V_Out;zeros(length(Z(:,1)),1)];
        V_Out >= Demand;
        V_In >=0; V_Out >=0; V>=0;
        V <= Branch(:,cap);
        ];
ops = sdpsettings('solver','gurobi','verbose',1);
solvesdp(Cons,Cost,ops);

result.V = double(V)
result.Cost = double(Cost)

    

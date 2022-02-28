%PSO >> function for the PSO ALGORITHM
%
% USAGES:   1.) [fxmin, xmin, Swarm, history] = PSO(psoOptions);
%           2.) [fxmin, xmin, Swarm, history] = PSO;
%           3.) fxmin = PSO(psoOptions);
%           3.) PSO 
%           etc.
%
% Arguments     : psoOptions--> A Matlab stucture containing all PSO related options. (see also: get_psoOptions)
% Return Values : [fxmin, xmin, Swarm, history]
%                     |     |       |       |_The history of the algorithm. Depicts how the function value of GBest changes over the run.
%                     |     |       |_The final Swarm. (A matrix containing co-ordinates of all particles)
%                     |     |_The co-ordinates of Best (ever) particle found during the PSO's run.
%                     |__The objective value of Best (^xmin) particle.
%
%  History        :   Author      :   JAG (Jagatpreet Singh)
%                     Created on  :   05022003 (Friday. 2nd May, 2003)
%                     Comments    :   The basic PSO algorithm.
%                     Modified on :   0710003 (Thursday. 10th July, 2003)
%                     Comments    :   It uses psoOptions structure now. More organized.
%  
%                     see also: get_psoOptions
function [fxmin, xmin, Swarm, accept_fes, success, history] = HMPSO(psoOptions,T,initSWARM,initStep)
%     need init_exemplers ,updating_exmplersCL ,updating_exmplersCLformerged,init_DLPE,update_DLPE
			%Initializations
			if nargin == 0
				psoOptions = get_psoOptions;
			end

			close all;
			%For Displaying 
			if psoOptions.Flags.ShowViz
				global vizAxes; %Use the specified axes if using GUI or create a new global if called from command window
				vizAxes = plot(0,0, '.');
				axis([-1000 1000 -1000 1000 -1000 1000]);   %Initially set to a cube of this size
				axis square;
				grid off;
				set(vizAxes,'EraseMode','xor','MarkerSize',15); %Set it to show particles.
				pause(1);
			end
			%End Display initialization
				 
			% 初始化种群
			Swarm = rand(psoOptions.Vars.SwarmSize, psoOptions.Vars.Dim)*(psoOptions.Obj.ub-psoOptions.Obj.lb) + psoOptions.Obj.lb; 
			me = psoOptions.Vars.max_gen;
			[SwarmSize,Dim]=size(Swarm);
		    % Inertia parameter : Weight C1 C2
			Weight = 0.9-(1:me)*0.5/me;
			 Weight1 = 1.5-(1:me)*1/me;
            CR = repmat(0.8,1,Dim);  %% 初始化重置概率
			c1 = 1.49445.*ones(me,1);
			c2 = 1.49445.*ones(me,1);
			flowrate = 5;
			topsize_ini = round(SwarmSize/flowrate);
			% CEC parameter
			[M1,M2,shifto,lambda10,lambda100]= functionparameter(SwarmSize,Dim);
			
            % 统计FES 
			success = 0; % Success Flag
			f2eval = psoOptions.Obj.f2eval;
			FEs = 0; % Function evaluations' counter
			G = 0; % 迭代次数
			k=0;
            accuracy = psoOptions.Vars.Threshold;
            maxruntime = psoOptions.Vars.Dim *10000; % 以运行次数计算为标准
            
			% init global
			fSwarm = feval(f2eval,Swarm,M1,M2,shifto,lambda10,lambda100);
			% fSwarm = feval(f2eval,Swarm);
			FEs = FEs + SwarmSize;
			pbest_pos = Swarm;
			pbest_val= fSwarm;
			[gbest_val,gbest_id] = min(pbest_val);
			gbest_pos = pbest_pos(gbest_id,:);
			history = [FEs, gbest_val];
			
			%% 划分子群, 以相对距离与适应度为标准设置标签
			[subSwarm,cl] = divswarm_fortest(Swarm,fSwarm,M1,M2,shifto,lambda10,lambda100, topsize_ini);
			subSwarmSize = size(subSwarm,2);
			%%%  设定速度边界
			VRmin = psoOptions.Obj.lb;
			VRmax = psoOptions.Obj.ub;
			mv=0.2*(VRmax-VRmin);
			v_min=repmat(-mv,SwarmSize,Dim);
			v_max=-v_min;
			for i=1:subSwarmSize
			 
				 pos_g{i} = subSwarm{i};
				 num_g(i) = size(pos_g{i},1);
				 % vel_g{i} = v_min(1:num_g(i),:) + (v_max(1:num_g(i),:) - v_min(1:num_g(i),:)).*rand(num_g(i),Dim);
				 vel_g{i}= v_min(1:num_g(i),:)+2.*v_max(1:num_g(i),:).*rand(num_g(i),Dim);    
				 % 初始化每个子群的最优值
				 pbest_pos_g{i} = pos_g{i};
				 pbest_val_g{i} = pbest_val(find(cl==i),:);
				 [gbest_val_g(i),gbest_id(i)] = min(pbest_val_g{i});
				 gbest_pos_g{i} = pbest_pos_g{i}(gbest_id(i),:);
				 gbest_val_g_r(i) = gbest_val_g(i);	
				 stop{i} = zeros(num_g(i),1); ;
				 %%% 储存每个子群需要重置的粒子信息 %%%
				 [gworst_val_g(i),gworst_id(i)] = max(pbest_val_g{i});
				 gworst_pos_g{i} = pbest_pos_g{i}(gworst_id(i),:);
				 %%% 用于计算每个维度的重置次数 %%%
                 CR_count{i}=zeros(num_g(i),Dim);
				 
				%%%%  邻居 %%%%%%%%
					 % 检查每个粒子是否更新
					 obj_func_slope_g{i}= zeros(num_g(i),1); 
					 num1=max(num_g);
					 t=0:1/(num_g(i) - 1):1;t=5.*t;
					 Pc{i}=0.0+(0.5-0.0).*(exp(t)-exp(t(1)))./(exp(t(num_g(i)))-exp(t(1)));
					 % 初始化每个粒子的最好邻居
					 fri_best_pos_g{i} = init_fri(  pbest_pos_g{i} ,pbest_val_g{i}',Pc{i},num_g(i),Dim);
				 
			
			end
		
		
		if psoOptions.Disp.Interval & (rem(k, psoOptions.Disp.Interval) == 0)
			fprintf('Iterations\t\tfGBest\t\t\tfevals\n');
		end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  THE  PSO  LOOP                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while FEs < maxruntime
         
        G = G + 1;
	    k = k + 1;	 
		if k > me
			 k = k -1;
		end
%%%%%%%%%%%%%%%%%%%%%%%  update subswarm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	    
	
	for i=1:subSwarmSize
	    num_g(i) = size(pos_g{i},1);
	    %%%%%%%%%%  更新邻居  %%%%%%%
		
		[fri_best_pos_g{i},obj_func_slope_g{i}] = updating_fri(fri_best_pos_g{i},pbest_pos_g{i} ,pbest_val_g{i}',Pc{i},num_g(i),Dim,obj_func_slope_g{i}); 
        A = repmat(gbest_pos, SwarmSize, 1);  %% 全局最优
		%%%%%%%%% different stategies to update velocity and position  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if i == 1
			vel_g{i} =  Weight(k) .* vel_g{i} + 1.49445 .* rand(num_g(i),Dim).*(fri_best_pos_g{i} - pos_g{i}) + 1.*1.49445.*rand(num_g(i),Dim).* (A(1:num_g(i),:)-pos_g{i});
%           vel_g{i} = 0.5 .* vel_g{i} + 1.49445 .*rand(num_g(i),Dim).*(pbest_pos_g{i}- pos_g{i}) + 1.49445 .*rand(num_g(i),Dim).* (A(1:num_g(i),:)-pos_g{i});
        else
			vel_g{i} = Weight(k) .* vel_g{i} + 1.49445 .* rand(num_g(i),Dim).*(fri_best_pos_g{i} - pos_g{i}); %邻域型
		end
		%检查边界
		vel_g{i}=(vel_g{i} < v_min(1:num_g(i),:)).*v_min(1:num_g(i),:)+(vel_g{i} >= v_min(1:num_g(i),:)).*vel_g{i};
		vel_g{i}= (vel_g{i} > v_max(1:num_g(i),:)).*v_max(1:num_g(i),:)+(vel_g{i} <= v_max(1:num_g(i),:)).*vel_g{i};
		%更新位置
		pos_g{i} = pos_g{i} + vel_g{i};	
		%更新局部与全局最优,fSwarm是列向量
		fSwarm_g{i} = feval(f2eval,pos_g{i},M1,M2,shifto,lambda10,lambda100);
		FEs = FEs + num_g(i); 
		% Updating the best position for each particle
		changeRows = [];
		changeRows = fSwarm_g{i} < pbest_val_g{i};
		pbest_val_g{i}(find(changeRows)) = fSwarm_g{i}(find(changeRows));
		pbest_pos_g{i}(find(changeRows), :) = pos_g{i}(find(changeRows), :);
		[gbest_val_g(i),gbest_id(i)] = min(pbest_val_g{i});
		%% 记录停滞信息 %%%%%%%%%%%%%
		
		obj_func_slope_g{i}(find(~changeRows))=obj_func_slope_g{i}(find(~changeRows))+1;
		stop{i}(find(~changeRows)) = stop{i}(find(~changeRows))+1;
		stop{i}(find(changeRows)) = 0; 
		
		%% 更新子群最优  %%%%%%%%%%%%%%%
		if  gbest_val_g_r(i) > gbest_val_g(i);	
		    gbest_val_g_r(i) = gbest_val_g(i);	
			gbest_pos_g{i} = pbest_pos_g{i}(gbest_id(i),:);
		end	  
		
		%% 更新全局最优  %%%%%%%%%%%%%%%
		
		if gbest_val > gbest_val_g_r(i)
		   gbest_val = gbest_val_g_r(i);
		   gbest_pos = gbest_pos_g{i};
		end 
		   
			
		%%%%%%%%%%%%%%%%%%%%% 交叉重置策略 %%%%%%%%%%%%%%%%%%%%%%%%	
        %% 将每层多代不更新的的粒子与全局最优交叉，并按照适应度分配到不同层级 %%%

        resetposid = [];
        [gworst_val_g(i),gworst_id(i)] = max(pbest_val_g{i});
		[value_rank,ordval] = sort(pbest_val_g{i},'descend');				
		if i == 1
			resetposid = ordval(1);  
		else
			resetposid = ordval(1);     %可调节每个子群重置粒子个数
        end
        
		for j = 1:size(resetposid)
            
            if  stop{i}(resetposid(j)) > 30
                %%% 1.和全局最优粒子随机交叉 %%%%
                gworst_pos_g{i} = pbest_pos_g{i}(resetposid(j),:);
                CR_used = CR - 0.05.*CR_count{i}(resetposid(j),:);
                changecr = rand(1,Dim) < CR_used;
                gworst_pos_g{i}(find(changecr)) = gbest_pos(find(changecr)); 
                %%%%%%%%%%%%%%%%  如果结果较好，则采用%%%%%%%%%%%%%
                changeval = feval(f2eval,gworst_pos_g{i},M1,M2,shifto,lambda10,lambda100);   
                FEs = FEs + 1; 
                if changeval < gworst_val_g(i)
                   CR_count{i}(resetposid(j),find(changecr))= CR_count{i}(resetposid(j),find(changecr))+1; 
                   pbest_pos_g{i}(resetposid(j),:) = gworst_pos_g{i};
                   %策略有效的话，则将重置后的粒子滞留值置零
                   stop{i}(resetposid(j)) = 0;
                end   
                %%% 若重置粒子优于全局最优且位于底层群，则将其与顶层粒子中最差的交换
                if i==2&&changeval < gbest_val
                   gbest_val = changeval;
                   gbest_pos = gworst_pos_g{i};
                   pbest_pos_g{i}(resetposid(j),:) = pbest_pos_g{1} (gworst_id(1),:);
                   pbest_pos_g{1} (gworst_id(1),:) = gworst_pos_g{i};
                end   
               %%%% 当该维度重置概率小于0.3时，将其重置次数归零
                   CR_count{i}(CR_count{i}>8)=0; 
            end
	     end		

		
	end
		%%%%%%%%%%%%%%%%%%%%% 流动策略 %%%%%%%%%%%%%%%%%%%%%%%%
		%% 将底层与中层粒子最好的粒子慢慢流向顶层，以提高后期的开采能力 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if FEs > maxruntime/2 && (mod(G, 1000)==0) && num_g(2) > topsize_ini
		clear swarm vel pbest_pos pbest_val;
		   Swarm = [pos_g{1};pos_g{2}];
		   vel = [vel_g{1};vel_g{2}];
		   pbest_pos = [pbest_pos_g{1};pbest_pos_g{2}];
		   pbest_val = [pbest_val_g{1};pbest_val_g{2}];
		
		   fprintf(' 粒子流动,底层群规模为 %d ',num_g(2));
		   %%% 将底层较优的粒子慢慢流向顶层%%%
			[subSwarm,cl] = divswarm_fortest(Swarm,pbest_val',M1,M2,shifto,lambda10,lambda100,num_g(1)+topsize_ini);

		    for i=1:subSwarmSize			 

				 pos_g{i} = subSwarm{i};
				 num_g(i) = size(pos_g{i},1);
				 % vel_g{i} = v_min(1:num_g(i),:) + (v_max(1:num_g(i),:) - v_min(1:num_g(i),:)).*rand(num_g(i),Dim);
				 lab = find(cl==i);
				 vel_g{i}= vel(lab,:);
				 pbest_pos_g{i} = pbest_pos(lab,:);
				 pbest_val_g{i} = pbest_val(lab,:);
				 [gbest_val_g(i),gbest_id(i)] = min(pbest_val_g{i}); 
				 stop{i} = zeros(num_g(i),1); ;
				 %%% 用于计算每个维度的重置次数 %%%
                 CR_count{i}=zeros(num_g(i),Dim);
				 
				%%%%  邻居 %%%%%%%%
					 % 检查每个粒子是否更新
					 obj_func_slope_g{i}= zeros(num_g(i),1); 
					 num1=max(num_g);
					 t=0:1/(num_g(i) - 1):1;t=5.*t;
					 Pc{i}=0.0+(0.5-0.0).*(exp(t)-exp(t(1)))./(exp(t(num_g(i)))-exp(t(1)));
					 % 初始化每个粒子的最好邻居
					 fri_best_pos_g{i} = init_fri(pbest_pos_g{i} ,pbest_val_g{i}',Pc{i},num_g(i),Dim);
			end
		   
		 
		end	
	       
   
    if success == 0
            if gbest_val <= accuracy
                accept_fes = FEs;
                success = 1;
            else
                accept_fes = 0;
            end
    end
        
    if psoOptions.Save.Interval & (rem(FEs, psoOptions.Save.Interval) == 0)
        history((size(history,1)+1), :) = [FEs, gbest_val];
    end
	
    if psoOptions.Disp.Interval & (rem(k, psoOptions.Disp.Interval) == 0)
        fprintf('%4d\t\t\t%.5g\t\t\t%5d\t\t\t%5d\n', G, gbest_val, FEs,T);
    end
	
end	
	
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  END  OF PSO  LOOP                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fxmin = gbest_val;
xmin = gbest_pos;
Swarm = pbest_pos;
  history = history(1:maxruntime/(3.*SwarmSize),:);

function [stdq]=StdQstandard(q,W);
while(find(abs(q)>W/2)>0)
    q(find(abs(q)>W/2))=abs(W/2-q(find(abs(q)>W/2)));
end
stdq=q;
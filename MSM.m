function Pc=MSM()  
%This function builds a Markov state model for two different example systems:
%model 1: macrophage polarization (https://doi.org/10.1016/j.jtbi.2023.111634) 
%model 2: bistable toggle switch
%Code written by Susanna Röblitz (susanna.roblitz@uib.no) and Maryam Yousefian (maryam.yousefian@uib.no)

    close(); 

    addpath("Dirichlet/")
    addpath("PCCA/")

    seed=2;
    rng(seed);
         
    model=input('Select a model (type 2 for toggle switch, type 1 for macrophage, ):');
    
    if model==1
        macrophage_exp3
        tf=10;
    elseif model==2
        tf=5000;
        toggle_exp1
    else
        error('invalid model number')
    end     

    %generate Voronoi cells and remove dublicated centers  
    ii=1;
    centers=[]; 
    while ii<= N
        %ii
        rx=round(Nx*rand(1));
        ry=round(Ny*rand(1));
        centers(ii,:) = [rx,ry];
        centers=unique(centers,'rows');
        ii=size(centers,1)+1;
     end


    N=size(centers,1)
    
    %create table with Voronoi cells to track their refinement later on
    index=[1:N]';
    log_ans= zeros(N,1);
    child1=zeros(N,1);
    child2=zeros(N,1);
    exit=zeros(N,1);
    voronoi_table = table(index,centers,log_ans,child1,child2,exit);
    
    

    %plot Voronoi tesselation
    figure(2);
    voronoi(centers(:,1),centers(:,2));
    xlim([0 Nx]);
    ylim([0 Ny]);
    title 'Voronoi tessellation diagram';
    hold on
    
    hor_sample_list=1:N;  %list with cells to be sampled horizontally
    level=1;    %set hierarchy level equal to 1
    centers=centers';
    while ~isempty(hor_sample_list) && level<=max_level
        refine_list_not_converge=[];
        %list_converge=[];
        mean_exit=zeros(1,length(hor_sample_list));
        for k=hor_sample_list
            k
            X=zeros(max_rep*N_ssa,D,num_traj);  %array holding all horizontal sampling points
            R=0;    %Gelman-Rubin convergence indicator
            Xtmp=[];
            rep=0;
            x0=voronoi_table.centers(k,:)';
            x0_vec=zeros(D,num_traj);
            for traj=1:num_traj
                x0_vec(:,traj)=x0;
            end
            t_exit_arr=zeros(1,num_traj);   %exit times for all rajectories
            %perform SSA within cell k (horizontal sampling):
            t=0;
            while (abs(R-1)>Rthr) && (rep<max_rep)
                rep=rep+1;
                for traj=1:num_traj   %generate several trajectories
                    %traj
                    if model==1
                        [Xt,tvec,t_exit]=macrophage_SSA(x0_vec(:,traj),N_ssa,0,k,voronoi_table,N);
                    elseif model==2
                        [Xt,tvec,t_exit]=toggle_SSA(x0_vec(:,traj),N_ssa,0,k,voronoi_table,N);
                    end
                    x0_vec(:,traj)=Xt(:,end);   %set lat state as new initial state
                    X((rep-1)*N_ssa+1:rep*N_ssa,:,traj)=Xt';
                    Xtmp=[Xtmp,Xt]; %concatinate al horizontal trajectories
                    %trajectory has left the cell for the first time
                    if ~isinf(t_exit) && t_exit_arr(traj)==0
                        t_exit_arr(traj)=t_exit+t;  %save first exit time
                    end
                    t=t+tvec(end);  %update time
                end
                %compute scale reduction factor based on 2nd half of traj
                start=ceil(rep*N_ssa/2);
                R=max(psrf(X(start:rep*N_ssa,:,:)));    
            end
            %compute mean exit time
            if min(t_exit_arr)==0
                 mean_exit(k)=Inf;
            else
                 mean_exit(k) = mean(t_exit_arr);
                 voronoi_table.exit(k)=mean(t_exit_arr);
            end
            if abs(R-1)>Rthr
                sprintf('convergence was not achieved for voronoi cell %d', k)
                refine_list_not_converge=[refine_list_not_converge,k];
%             else
%                 sprintf('convergence was achieved for voronoi cell %d', k)
%                 list_converge=[list_converge,k];
            end
            field_is = strcat('k', num2str(k));
            Xh.(field_is)=Xtmp;
%             if level>1
%                  figure(2)
%                  plot(Xtmp(1,:),Xtmp(2,:),'*')  
%                  hold on
%             end
         
        end %for k=hor_sample_list

        %create list of cells to be refined
        [~,idx]=sort(mean_exit(refine_list_not_converge),'descend');
        final_refine_list= refine_list_not_converge(idx);
        final_refine_list=final_refine_list(1:min(max_ref,length(final_refine_list)));
        %add new Voronoi cells
        if level<max_level
            voronoi_table.log_ans(final_refine_list)=logical(1);
            %extend the table of nodes
            n_new=2*length(final_refine_list);  %length of new table
            index=[N_total+1:N_total+n_new]';
            log_ans=logical(zeros(n_new,1));
            centers=zeros(n_new,2);
            child1=zeros(n_new,1);
            child2=zeros(n_new,1);
            exit=zeros(n_new,1);
            T_new=table(index,centers,log_ans,child1,child2,exit);
            voronoi_table=[voronoi_table;T_new];

            counter= N_total;
            %create two new children nodes for every cell that is refined
            for ff=final_refine_list
                voronoi_table.child1(ff)=counter+1;
                voronoi_table.child2(ff)=counter+2;
                fieldname = strcat('k', num2str(ff));
                XR = Xh.(fieldname)';   
                [idx,C1] = kmeans(XR,2);
                voronoi_table.centers(counter+1,:)=round(C1(1,:));
                voronoi_table.centers(counter+2,:)=round(C1(2,:));
                counter=counter+2;
                figure(2);
                plot(C1(:,1),C1(:,2) ,'kx','MarkerSize',5,'LineWidth',1)
                %title 'Cluster Assignments and Centroids'
                hold on
             end    %for ff=final_refine_list  
             % add new cells to list for horizontal sampling
             hor_sample_list=(N_total+1:N_total+n_new);
             N_total=N_total+n_new; %new total number of cells    
        end  %if level<max_level
        level=level+1
    end    %end horizontal sampling

    

    %save results from horizontal sampling
    save_name = sprintf('model%d_seed%d_horSampling__N%d__N_ssa%d__N_vert0%d_tf%d.mat', model,seed,N, N_ssa, N_vert,tf);
    save(save_name,'voronoi_table','Xh')

    %determine tf for vertical sampling:
    %tmf=max(mean_exit);
    tmf=max(voronoi_table.exit(:))
    sprintf('maximum mean exit time: %f',tmf)

    sprintf('start vertical sampling:');
    cells=find(voronoi_table.log_ans==0);
    num_cells=length(cells);
    num_hor=0;  %counter for total number of horizontal sampling points
    counter=0;  %counter for cells 
    P=zeros(num_cells,num_cells);
    Pmle=zeros(num_cells,num_cells);
    alphaMat=zeros(num_cells,num_cells);
    for ss = cells'
        ss
        counter=counter+1;
        rows=zeros(N_rows,N_total);
        X_hor=Xh.(strcat('k', num2str(ss)));
        len_hor=size(X_hor,2);
        num_hor=num_hor+len_hor;
        for k=1:N_rows
            Xtmp=X_hor(:,randsample(1:len_hor,N_vert)); %sample N_vert points from the horizontal traj
            for jj=1:size(Xtmp,2) 
                if model==1
                    [Xt,~,~]=macrophage_SSA(Xtmp(:,jj),0,tf,0,voronoi_table,N);%height(voronoi_table));
                elseif model==2
                    [Xt,~]=toggle_SSA(Xtmp(:,jj),0,tf,0,voronoi_table,N);%height(voronoi_table));
                end
                cell_idx = membership(Xt(:,end),voronoi_table,N);   %check in which cell the traj ends
                rows(k,cell_idx)=rows(k,cell_idx)+1;                 
            end
        end %for k=1:N_rows
        rows=rows(:,cells);
        Pmle(counter,:)=sum(rows);  %maximum likelihood estimate
        rows=diag(1./sum(rows,2))*rows; %normalize rows
        alpha=dirichlet_fit_newton(rows);
        sig=sum(alpha);
        mu=alpha/sig;
        var=mu.*(1-mu)/(sig+1);
        alphaMat(counter,:)=alpha;
        P(counter,:)=mu;    %Dirichlet estimate
    end %for ss=cells'  
    %Pmle=diag(1./sum(Pmle,2))*Pmle;  

    %num_total_points=num_hor+N_total*N_vert;

    save_name = sprintf('model%d_seed%d_alphaMat__N%d__N_ssa%d__N_vert0%d_tf%d.mat', model,seed,N, N_ssa, N_vert,tf);
    save(save_name,'alphaMat','Pmle')

    stop=1;

      

end

%this file produces the figures for the toggle switch model (model A)

addpath("PCCA/")
addpath("Dirichlet/")

seed=2;
rng(seed);

Nx=300; Ny=300;
nc=2;
%toggle_refsol();

%case A1
    load('Results/model2_seed2_alphaMat__N100__N_ssa500__N_vert0100.mat');
    load('Results/model2_seed2_horSampling__N100__N_ssa500__N_vert0100.mat');
%case A2    
    %load('Results/model2_seed2_alphaMat__N100__N_ssa500__N_vert0500.mat');
    %load('Results/model2_seed2_horSampling__N100__N_ssa500__N_vert0500.mat');
%case A3
    %load('Results/model2_seed2_alphaMat__N20__N_ssa500__N_vert0500.mat');
    %load('Results/model2_seed2_horSampling__N20__N_ssa500__N_vert0500.mat');

    
    target=1+eps;               %target eigenvalue

    cells=find(voronoi_table.log_ans==0);
    cell_centers=voronoi_table.centers(cells,:);
    
    %compute average MSM model (maximum likelihood)
    Pmle=diag(1./sum(Pmle,2))*Pmle;
    [Pc,chi,A,ovec]=computeMSM(Pmle,nc,target,cell_centers',Nx,Ny);

    %compute MSM with the mean of the Dirichlet distribution
    %P=computeP(alphaMat);
    %[Pc,chi,A,ovec]=computeMSM(P,nc,target,cell_centers',Nx,Ny);

    %make chi coarse
    chic=zeros(size(chi,1),size(chi,2));
    [~,idx]=max(chi,[],2);
    for i=1:size(chic,1)
        chic(i,idx(i))=1;
    end 

    numMSM=20;
    %use chic as reference for sorting chi correctly
    [Pmats,chimats,eigvals,weights]=computeMSM_with_error_dirichlet(alphaMat,nc,target,numMSM,cell_centers',Nx,Ny,chic);
    %generate boxplots
    visualize_errors2(Pmats,chimats,eigvals,weights,numMSM,cell_centers',chic,Nx,Ny);  
    
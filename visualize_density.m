function visualize_density(rho,voronoi_table,N,points,Nx,Ny)
%visualize density of horizontal sampling points weighted with rho

    cells=find(voronoi_table.log_ans==0);
    cell_centers=voronoi_table.centers(cells,:);


    
    rho_matrix=zeros(Nx+1,Ny+1); %value of density in every grid point
    for i=0:Nx
        for j=0:Ny
            cell=membership([i;j],voronoi_table,N);
            idx=find(cells==cell);
            rho_matrix(i+1,j+1)=rho(idx);
        end
    end

   

    count_matrix=sparse(Nx+1,Ny+1);
    out_count=0;
    for ss = cells'
        ss
        X_hor=points.(strcat('k', num2str(ss)));
        num_hor=size(X_hor,2);
        for k=1:num_hor
            i=X_hor(1,k);
            j=X_hor(2,k);
            if i<(Nx+1) && j<(Ny+1)
                count_matrix(i+1,j+1)=count_matrix(i+1,j+1)+rho_matrix(i+1,j+1)/num_hor;
            else
                out_count=out_count+1;
            end
        end
    end

    %contour(Z): The column and row indices of Z are the x and y coordinates in the plane, respectively.
    contour(count_matrix')
    set(gca,'FontSize',18)
    axis equal
    xlim([0,Nx])
    xlim([0,Ny])
    colorbar()

    stop=1;

end

function [cell_idx] =membership_new(x,voronoi_table,N)
    
    centers=voronoi_table.centers(1:N,:);
    centers=centers';

    sc=(centers-x).^2;  
    edist=sc(1,:)+sc(2,:); 
    [~,cell_idx]=min(edist); 

    while voronoi_table.log_ans(cell_idx)==1
        child1=voronoi_table.child1(cell_idx);
        child2=voronoi_table.child2(cell_idx);
        children=[child1,child2];
        center1=voronoi_table.centers(child1,:);
        center2=voronoi_table.centers(child2,:);
        centers=[center1;center2]';
        sc=(centers-x).^2;  
        edist=sc(1,:)+sc(2,:); 
        [~,cell_idx]=min(edist);
        cell_idx=children(cell_idx);
    end
    
end
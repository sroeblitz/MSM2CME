function chiR=macrophage_Voronoi_chi(centers,Nx,Ny,chiref)
%compute reference membership vectors at center points

    nrows=size(centers,2);
    
    %load('macrophage_chiref.mat');  
    nc=size(chiref,2);
    %find values of membership functions at centerpoints
    %generate surface plots of chi
    %Ncells=size(centers,2);
    chiR=zeros(nrows,nc);
    for j=1:nrows
        %find value of chiref at centerpoint
        point=centers(:,j);
        stateIndex=point(1)*(Nx+1)+point(2)+1;
        chivalue=chiref(stateIndex,:);
        [~,column]=max(chivalue);
        chiR(j,column)=1;           
    end

%         %generate surface plots of chi
%         [Xq,Yq] = meshgrid(0:1:Nx,0:1:Ny);
%         for i=1:nc
%             figure(2+i)
%             %plot3(centers(1,:),centers(2,:),chi(:,i),'ro')
%             %hold on
%             F=scatteredInterpolant(centers(1,:)',centers(2,:)',chiR(:,i));
%             F.Method = 'natural';
%             vq=F(Xq,Yq);
%             %mesh(Xq,Yq,vq)
%             contour(Xq,Yq,vq)
%             title('chi',i)
%         end

end
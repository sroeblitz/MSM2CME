function idx=reorder_chic(chi,chiref)
%reorder the coarse chi so that it matches with chiref

    [n,nc]=size(chi);
    idx=zeros(1,nc);
    for i=1:nc
        [~,index]=min(sum((chi-chiref(:,i)).*(chi-chiref(:,i))));
        idx(i)=index;
    end
    if length(unique(idx))<nc
        warning('reordering of membership vectors went wrong')
    end
end


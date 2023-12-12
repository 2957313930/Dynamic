function [atf_data] = constraint(pre_data,bound,m,N)
    atf_data = zeros(size(pre_data));
    for i=1:1:m
        for_index = pre_data(i,:);
        indL = for_index<bound(i,1); % 超过约束下限的索引
        indU = for_index>bound(i,2); % 超过约束上限的索引
        % 对超过约束部分的值重新赋值
        a = bound(i,1) + (bound(i,2)-bound(i,1))*rand(1,N);
        b = bound(i,1) + (bound(i,2)-bound(i,1))*rand(1,N);
        for_index(indL) = a(indL);
        for_index(indU) = b(indU);
        atf_data(i,:) = for_index; % 完成约束赋值
    end
end

function [ tf ] = test_functions( pts )
%TEST_FUNCTIONS Calculate the mean of an assortment of test functions using
%particle filter approximations to the joint distributions

tf = cell(10,1);

for ii = 1:length(pts)
    
    tf{1} = sum(pts{ii}(1:2,:).^2, 1);
    tf{2} = sum(pts{ii}(3:4,:).^2, 1);
    tf{3} = sum(pts{ii}(1:2,1:end-1).*pts{ii}(1:2,2:end), 1);
    tf{4} = sum(pts{ii}(3:4,1:end-1).*pts{ii}(3:4,2:end), 1);
    
end

end


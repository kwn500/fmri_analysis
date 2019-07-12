function A = GetConvolutionMatrix(b,numA)
    A = zeros(numA,numA);
    vec = [b  zeros(1,numA-numel(b))];
    for i=1:size(A,1)
        A(i,:) = circshift(vec,[1 i]);
    end
end
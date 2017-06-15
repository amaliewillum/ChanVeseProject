function p = init_phi(im,type)
    im = dimensionz(im);
    %m = im(150:250,150:250);
    [dim1, dim2] = size(im);
    p = zeros(dim1+2,dim2+2);
    
    switch lower (type)
        case 'circle'
            for i = 1:dim1+2
                for j = 1:dim2+2
                    p(i,j) = (sqrt((i/dim1-0.5)^2 + (j/dim2-0.5)^2) - 0.2) * 30;
                end
            end
            
        case 'grid'
            for i = 1:dim1+1
                for j = 1:dim2+1
                    p(i,j) = sin(i*pi/5) + sin(j*pi/5);
                end
            end
            
        case 'circle 2'
            for i = 1:dim1+2
                for j = 1:dim2+2
                    p(i,j) = (sqrt(((i+2)/dim1-0.3)^2 + (j/dim2-0.9)^2) - 0.2) * 30;
                end
            end
            
        case 'square'
            %p = zeros(dim1+2,dim2+2);
            p(floor((dim1+2)/3:(dim1+2)*2/3),floor((dim2+2)/3:(dim2+2)*2/3)) = 1;  
%             p(floor(((dim1+2)/3)+1):floor((dim1+2)*2/3)-1,floor((dim2+2)/3)+1:floor((dim2+2)*2/3)-1) = 0;
            p = bwdist(p)-bwdist(1-p)+im2double(p)-.5;
    end
end
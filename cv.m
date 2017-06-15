function [] = cv(im,lambda1,lambda2,num_iter,mu,nu,dt,bw,j,phi0,Noize,DJ)
    
    format long  
    eps = 0.001;
    
    h = 0.1;
    epsilon = 1;
    
    P = dimensionz(im);

    I = P;
    image = P;
    image = image ./ max(image(:)+eps);
    image = image ./ mean2(image);
    
    
    
    switch Noize
        case 'Gauss'
            image = imnoise(image,'gaussian');
            I = image;
        case 'Speckle'
            image = imnoise(image,'speckle');
            I = image;
        case 'Salt & Pepper'
            image = imnoise(image,'salt & pepper',0.1);
            I = image;
        case 'Blur 1'
            H = fspecial('gaussian',size(image),1);
            image = imfilter(image,H,'replicate');
            I = image;
        case 'Blur 2'
            H = fspecial('gaussian',size(image),2);
            image = imfilter(image,H,'replicate');
            I = image;
        case 'Blur 3'
            H = fspecial('gaussian',size(image),3);
            image = imfilter(image,H,'replicate');
            I = image;
        case 'Blur 5'
            H = fspecial('gaussian',size(image),5);
            image = imfilter(image,H,'replicate');
            I = image;
        case 'Blur 10'
            H = fspecial('gaussian',[30 30],10);
            image = imfilter(image,H,'replicate');
            I = image;
        case 'Blur 15'
            H = fspecial('disk',100);
            image = imfilter(image,H,'replicate');
            I = image;
        case 'Blur + Gauss'
            H = fspecial('gaussian',size(image),5);
            image = imfilter(image,H,'replicate');
            image = imnoise(image,'speckle',1);
%             image = imnoise(image,'gaussian');
            I = image;
        otherwise
    end
  
    phi = phi0;
    phi = phi./max(abs(phi(:))+eps);

    H = @(z,epsilon) 0.5.*(1+(2./pi).*atan(z./epsilon));
    delta = @(z,epsilon) 1/pi .* epsilon./(epsilon^2 + z.^2);


    MM = @(k) k(2:end-1,2:end-1);   
    Mp01 = @(k) k(1:end-2,2:end-1); 
    Mm10 = @(k) k(2:end-1,1:end-2); 
    Mm01 = @(k) k(3:end,2:end-1);   
    Mp10 = @(k) k(2:end-1,3:end);
    Mmm = @(k) k(3:end,1:end-2);
    Mpp = @(k) k(1:end-2,3:end);
    Mmp = @(k) k(1:end-2,1:end-2);
    Mpm = @(k) k(3:end,3:end);
    
    if bw
        c1 = @(p,epsilon) sum(image(:).*H(p(:),epsilon))./(sum(H(p(:),epsilon)));
        c2 = @(p,epsilon) sum(image(:).*(1-H(p(:),epsilon)))./(sum(1-H(p(:),epsilon)));
    else
        c1 = @(p,epsilon) 0;
        c2 = @(p,epsilon) 1;
    end
   
    time = 0;
    [dim1, dim2] = size(phi);
    figure('units','centimeters','position',[10 10 35 35]);
    margins = [.05 .05];
    subplot_tight(2,2,1,margins),imagesc(I),axis image, axis off, title('Original Image'), colormap(gray);
    hold on
    subplot_tight(2,2,2,margins),contour(flipud(phi0),[0 0], 'r'),axis image, axis off,  title('Initial Contour');
    
    subplot_tight(2,2,3,margins), showphi(image,MM(phi),0,P,time);
    subplot_tight(2,2,4,margins), seg(MM(phi),epsilon,c1(MM(phi),epsilon),c2(MM(phi),epsilon));
    
    
    phi_check = zeros(dim1,dim2,2);
    
        
    for i = 1:num_iter
        
        tic;
        
        Dxx = (Mp10(phi) - 2.*MM(phi) + Mm10(phi))./(h^2+eps);
        Dyy = (Mp01(phi) - 2.*MM(phi) + Mm01(phi))./(h^2+eps);
        Dxy = (Mpp(phi) - Mpm(phi) - Mmp(phi) + Mmm(phi))./(eps+4.*h.^2);%--++
        Dx = (Mp10(phi) - Mm10(phi))./(2.*h+eps);
        Dy = (Mp01(phi) - Mm01(phi))./(2.*h+eps); %-
        g = sqrt((Dx).^2 + (Dy).^2);

        kappa = (Dx.^2 .* Dyy + Dy.^2 .* Dxx - 2.*Dxy .* Dx .* Dy) ./ (eps + g.^3);
        kappa = kappa./max(abs(kappa(:))+eps);
        

        old_phi = phi;
        old_phi = old_phi./max(abs(old_phi(:))+eps);
        
        if mod(i,43) == 0
            phi_temp = phi_check(:,:,2);
            phi_check(:,:,2) = phi;
            phi_check(:,:,1) = phi_temp;
        end
        
        
        F = (-lambda1.*(image - c1(MM(phi),epsilon)).^2 ... 
            +lambda2.*(image - c2(MM(phi),epsilon)).^2); %./max(abs(F(:)))
        
        phi(2:end-1,2:end-1) = (F./(max(abs(F(:))) +eps) + mu.*kappa + nu).*delta(MM(phi),epsilon); 

           
        phi = phi ./ max(abs(phi(:))+eps);
        phi = old_phi + dt * phi;
        phi(1,1:end) = phi(2,1:end);
        phi(end,1:end) = phi(end-1,1:end);
        phi(1:end,1) = phi(1:end,2);
        phi(1:end,end) = phi(1:end,end-1);        
                
        time = time + toc;
        
%         phi_STOP = phi_check(:,:,1);
%         if isequal(phi > 0,phi_STOP > 0)
%             subplot_tight(2,2,3,margins), showphi(image,MM(phi),i,P,time);
%             subplot_tight(2,2,4,margins), seg(MM(phi),epsilon,c1(MM(phi),epsilon),c2(MM(phi),epsilon));
%             break; 
%         end
        
        if(mod(i,j) == 0)
            subplot_tight(2,2,3,margins), showphi(image,MM(phi),i,P,time);
            subplot_tight(2,2,4,margins), seg(MM(phi),epsilon,c1(MM(phi),epsilon),c2(MM(phi),epsilon));
        end
    
    end
   
    hold off

    if DJ
        inside_P = find(P == 255);
        if c1(MM(phi),epsilon) > c2(MM(phi),epsilon)
%             inside_phi = find(MM(phi) <= 0);
            inside_phi = find(MM(phi) >= 0);
        else
            inside_phi = find(MM(phi) <= 0);
%             inside_phi = find(MM(phi) >= 0);
        end

        QS = (2 * length(intersect(inside_P,inside_phi))) / (length(inside_P) + length(inside_phi));
        disp(['QS = ' num2str(QS)])

        J = length(intersect(inside_P,inside_phi)) / (length(union(inside_P,inside_phi)));
        disp(['J(A,B) = ' num2str(J)])
    end
end

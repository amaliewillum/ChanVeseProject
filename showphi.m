function phi_img = showphi(u0, phi0, i,p,Time)
    
    %imshow(u0,'InitialMagnification','fit');
%     u0 = u0 * max(p(:)-min(p(:))+eps) + min(p(:));
    u0 = u0 * max(p(:)) * mean(p(:));
    imagesc(u0), axis image, axis off, colormap(gray);
    hold on;
    contour(phi0, [0 0],'m','LineWidth',5);
    contour(phi0, [0 0],'c','LineWidth',2);
    hold off;
    title([num2str(i) ' Iterations and ' num2str(round(Time,2)) ' s CPU time']);
    drawnow;
end


function segmentation = seg(phi0,eps,C1,C2)
    imagesc((phi0 > 0).*C1 + (phi0 < 0).*C2), axis image, axis off,title('Region Based Segmentation');
%     imagesc(phi0 < 0), axis image, axis off,title('Region Based Segmentation');
    drawnow;
end
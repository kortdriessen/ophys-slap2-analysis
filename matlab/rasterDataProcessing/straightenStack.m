function [IMout, IMc, IMsk]= straightenStack(IM, IMc, IMsk)
Y = makeHighPass(sum(IM,3)); % data order is XYCZ
Z2fft = fft2(Y(:,:,1,1));
IMout = IM;
dX = 0; dY = 0;
for Z1 = 1:(size(IM, 4)-1)
    Z1fft = Z2fft;
    Z2fft = fft2(Y(:,:,1,Z1+1));
    output = dftregistration_clipped(Z1fft,Z2fft,4,8);
    dX = dX + output(4); dY = dY+output(3);
    IMout(:,:,:,Z1+1) = imtranslate(IM(:,:,:,Z1+1), [dX dY]);
    IMc(:,:,:,Z1+1) = imtranslate(IMc(:,:,:,Z1+1), [dX dY]);
    IMsk(:,:,:,Z1+1) = imtranslate(IMsk(:,:,:,Z1+1), [dX dY]);
end
end
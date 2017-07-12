function [ I ] = calc_int( v,n,img, rng )
%CALC_INT calculates the intensity profiles along the normal vectors (n) at
%points mp in the image (img).

% Make array of x-, y- and z-coordinates of points along normal vectors
X = v(:,1)*ones(1,length(rng)) + n(:,1) * rng;
Y = v(:,2)*ones(1,length(rng)) + n(:,2) * rng;

switch ndims(img)
    case 2
        % One channel in the data.
        I = interp2(double(img),X,Y);
    case 3
        % Multiple channels in the data. Interpolate each channel
        % independently.
        I = zeros(size(X,1),size(X,2),size(img,3));
        for ii = 1 : size(img,3)
            I(:,:,ii) = interp2(double(img(:,:,ii)),X,Y);
        end
end

end


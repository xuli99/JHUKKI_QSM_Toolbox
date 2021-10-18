function res = adjD(y)

%y1 = ones(imsize)*y(1)/sqrt(prod(imsize));
%yx = (reshape(y(2:prod(imsize)+1), imsize(1), imsize(2)));
%yy = (reshape(y(prod(imsize)+2:end), imsize(1), imsize(2)));

res = adjDx(y(:,:,:,1)) + adjDy(y(:,:,:,2)) + adjDz(y(:,:,:,3));

return;


function res = adjDz(x)
res = x(:,:,[1,1:end-1]) - x;
res(:,:,1) = -x(:,:,1);
res(:,:,end) = x(:,:,end-1);

function res = adjDy(x)
res = x(:,[1,1:end-1],:) - x;
res(:,1,:) = -x(:,1,:);
res(:,end,:) = x(:,end-1,:);

function res = adjDx(x)
res = x([1,1:end-1],:,:) - x;
res(1,:,:) = -x(1,:,:);
res(end,:,:) = x(end-1,:,:);



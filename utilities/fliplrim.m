function im = fliplrim(im)
im(:,:,1)=fliplr(im(:,:,1)); 
im(:,:,2)=fliplr(im(:,:,2)); 
im(:,:,3)=fliplr(im(:,:,3)); 
end
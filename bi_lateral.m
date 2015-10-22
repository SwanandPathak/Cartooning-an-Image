clear
close all
imtool close all
clc

img=(imread('Castle.jpg'));                                               %Reading image
startImg=img;
img=imresize(img,0.2);                                                   %resizing to 20% of original image. Do only for very large images to lessen computations.
img=im2double(img);                                                      %Converting to double data type
r=img(:,:,1);
g=img(:,:,2);                                                            %Separating R,G,B channels
b=img(:,:,3);
filt_size=13;                                                            %Setting filter size to 13
sigma_gau=3;                                                             %Setting sigma for normal gaussian to 3 
sigma_r=0.08;                                                            %Setting sigma for bilateral gaussian to 0.08 
gau=fspecial('gaussian',filt_size,sigma_gau);                            %Computing gaussian Filter, for given size and sigma
[rows cols ch]=size(img);                                                     
start=((filt_size-1)/2)+1;
new_img=img;
for(i=start:(rows-start+1))
    for(j=start:(cols-start+1))
        block=img((i-(start-1):i+(start-1)),(j-(start-1):j+(start-1)),:);  %Extracting block from image, depending on image size, having center pixel co-ords= (i,j)
        dif_r=(block(:,:,1)-img(i,j,1)).^2;
        dif_g=(block(:,:,2)-img(i,j,2)).^2;                                %Computing the Square of the difference values with respect to current center pixel for all channels
        dif_b=(block(:,:,3)-img(i,j,3)).^2;
        loc_gau_r=exp(-1*dif_r/(2*sigma_r*sigma_r));
        loc_gau_g=exp(-1*dif_g/(2*sigma_r*sigma_r));                       %Implemention the Bilateral Gaussian function for all channels
        loc_gau_b=exp(-1*dif_b/(2*sigma_r*sigma_r)); 
        prod(:,:,1)=block(:,:,1).*gau.*loc_gau_r;
        prod(:,:,2)=block(:,:,2).*gau.*loc_gau_g;                          %Element by element multiplication of all three-normal gauss,bilateral gauss,block
        prod(:,:,3)=block(:,:,3).*gau.*loc_gau_b; 
        new_img(i,j,1)=sum(sum(prod(:,:,1)));
        new_img(i,j,2)=sum(sum(prod(:,:,2)));                              %Combinating all 3 channels
        new_img(i,j,3)=sum(sum(prod(:,:,3)));           
    end
end
try
    img=rgb2gray(img);                                                     %converting the image to gray scale
end
img=im2double(img);                                                        %converting to double data type
h = fspecial('gaussian', 5, 1.4);                                          %Applying the gaussian filter to remove noise
w = imfilter(img,h,'same');
kgx=[-1,0,1;-2,0,2;-1,0,1];                                                %sobel-kernel
kgy=[1,2,1;0,0,0;-1,-2,-1];
wx=imfilter(w,kgx,'same');                                                 %directional gradients
wy=imfilter(w,kgy,'same');
[m,n]=size(wx);
g=zeros(m,n);
theta=zeros(m,n);
thetaR=zeros(m,n);
g=sqrt(wx.*wx+wy.*wy);                                                     %calculating gradient magnitude
theta=atand(wy./wx);                                                       %calculating tan-arc to get gradient direction
for i=1:m
    for j=1:n
       
        if(theta(i,j)>-22.5 && theta(i,j)<22.5)
            thetaR(i,j)=0;
        elseif(theta(i,j)>22.5 && theta(i,j)<67.5)                         %separating directions
            thetaR(i,j)=45;
        elseif(theta(i,j)>67.5 || theta(i,j)<-67.5)
            thetaR(i,j)=90;
        else
            thetaR(i,j)=135;
        end
    end
end
h=imhist(img,256);
h=h./(m*n);
cdf=cumsum(h);                                                             %Computing CDF
tstrong=(find(cdf>0.8,1,'first'))/255                                      %Selecting threshold which is the intensity corresponding to 80% value of the CDF
tweak=tstrong*0.4;                                                         %Weak threshold is 40% of Strong thresold
r=[];
c=[];
for i=1:m
    for j=1:n
        if(g(i,j)>tstrong)
            temp_s(i,j)=1;                  
            temp_w(i,j)=1;
            r=[r;i];
            c=[c;j];
        elseif(g(i,j)<tweak)
            temp_s(i,j)=0;
            temp_w(i,j)=0;
        else
            temp_w(i,j)=1;
            temp_s(i,j) = 0;
        end
       
    end
end
b=bwselect(temp_w,r,c,8);                                                  %Performing hysteresis 
b1=bwmorph(b,'thin','inf');                                                %Morphological erosion to get thin output
k=new_img;
for i=1:m
    for j=1:n
        if(b1(i,j)==1)
            k(i,j,1)=0;
            k(i,j,2)=0;                                                    %Adding the canny detected edge to the Bilateral filter output
            k(i,j,3)=0;
        end
    end
end
figure;
subplot(1,2,1)
imshow(startImg);
title('ORIGINAL INPUT IMAGE')
subplot(1,2,2)
imshow(k);
title('FINAL CARTOON IMAGE');

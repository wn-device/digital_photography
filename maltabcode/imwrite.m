path = "D:\work\data\small_hl.png";
write_path = "D:\work\data\small_hl.dat";
read_path = write_path;
close all;
img = imread(path);
% figure,imshow(img);
height =size(img,1);
width = size(img,2);
channel = size(img,3);
fp=fopen(write_path,"wb");
fwrite(fp,width,"uint32");
fwrite(fp,height,"uint32");
fwrite(fp,channel,"uint32");
for c=1:channel
    for h=1:height
       fwrite(fp,img(h,:,c),"float");
    end
end
fclose(fp);
fp=fopen(write_path,"rb");
width = fread(fp,1,"uint32");
height = fread(fp,1,"uint32");
channel = fread(fp,1,"uint32");
img_read_data = zeros(height,width,channel);
for c=1:channel
    for h=1:height
       img_read_data(h,:,c) = fread(fp,width,"float");
    end
end
fclose(fp);
figure,imshow(uint8(img_read_data));




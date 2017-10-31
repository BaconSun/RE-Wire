file_path='test\\';
filename='view';
fileidx='1';
A=imread([file_path filename fileidx '_bd0' '.jpg']);
B = ones(size(A,1), size(A,2));
thre=30;%eliminate the image error
sharp=2;%sharpness
for i=1:size(A,1)
    for j=1:size(A,2)
        if(A(i,j)>thre) 
            B(i,j)=0;
        end
    end
end

%Choice1:for real data;
im_ori=B;
im=B;
%im = bwmorph(im_ori,'close');
% im = imdilate(im_ori,strel('square',10));
%imshow(im);

% im = imerode(im,strel('square',sharp));
imshow(im);

%Choice2:for synthetic data
%im=B;

[edgelist, labelededgeim] = edgelink(im, 10);

drawedgelist(edgelist, size(im), 1, 'rand', 5); axis off
seglist = lineseg(edgelist, 2);
drawedgelist(seglist, size(im), 2, 'rand', 6); axis off
minlength=100;
nedgelist = cleanedgelist(edgelist, minlength);
drawedgelist(nedgelist, size(im), 1, 'rand', 7); axis off

fileID = fopen([file_path 'curvepoints' fileidx '.txt'], 'w');
fprintf(fileID, '%d\n', size(nedgelist, 2));


%output the thinning image ()
color=hsv(size(nedgelist,2));

color = color(randperm(size(nedgelist,2)),:)
segimg=zeros(size(A,1),size(A,2),3);
segimg2=zeros(size(A,1),size(A,2),3);
segimg(:,:,:)=1;
segimg2(:,:,:)=1;
for i=1:size(nedgelist,2)
for j=1:size(nedgelist{1,i})
segimg(nedgelist{1,i}(j,1),nedgelist{1,i}(j,2),:)=color(i,:);
segimg(nedgelist{1,i}(j,1)-1,nedgelist{1,i}(j,2)-1,:)=color(i,:);
segimg(nedgelist{1,i}(j,1)-1,nedgelist{1,i}(j,2),:)=color(i,:);
segimg(nedgelist{1,i}(j,1),nedgelist{1,i}(j,2)-1,:)=color(i,:);
segimg(nedgelist{1,i}(j,1)+1,nedgelist{1,i}(j,2),:)=color(i,:);
segimg(nedgelist{1,i}(j,1),nedgelist{1,i}(j,2)+1,:)=color(i,:);
segimg(nedgelist{1,i}(j,1)+1,nedgelist{1,i}(j,2)+1,:)=color(i,:);
segimg(nedgelist{1,i}(j,1)+1,nedgelist{1,i}(j,2)-1,:)=color(i,:);
segimg(nedgelist{1,i}(j,1)-1,nedgelist{1,i}(j,2)+1,:)=color(i,:);

segimg2(nedgelist{1,i}(j,1),nedgelist{1,i}(j,2),:)=0;
segimg2(nedgelist{1,i}(j,1)-1,nedgelist{1,i}(j,2)-1,:)=0;
segimg2(nedgelist{1,i}(j,1)-1,nedgelist{1,i}(j,2),:)=0;
segimg2(nedgelist{1,i}(j,1),nedgelist{1,i}(j,2)-1,:)=0;
segimg2(nedgelist{1,i}(j,1)+1,nedgelist{1,i}(j,2),:)=0;
segimg2(nedgelist{1,i}(j,1),nedgelist{1,i}(j,2)+1,:)=0;
segimg2(nedgelist{1,i}(j,1)+1,nedgelist{1,i}(j,2)+1,:)=0;
segimg2(nedgelist{1,i}(j,1)+1,nedgelist{1,i}(j,2)-1,:)=0;
segimg2(nedgelist{1,i}(j,1)-1,nedgelist{1,i}(j,2)+1,:)=0;

end
end
imshow(segimg);
imwrite(segimg, [file_path filename fileidx '_bd2' '.png']);
imwrite(segimg2, [file_path filename fileidx '_bd1' '.png']);



%colors=[0,0,1; 0,1,0; 1,0,0; 
%  1,1,0; 1,0,1; 0,1,1;
%  0,0,0.5; 0,0.5,0; 0.5,0,0; 
%  0.5,1,0; 1,0,0.5; 0,0.5,1;
%  1,0.5,0; 0.5,0,1; 0,1,0.5;
%  0,0.5,0.5; 0.5,0.5,0; 0.5,0,0.5;
%  1, 0.5, 0.5; 0.5, 1, 0.5; 0.5, 0.5, 1;];

%C = zeros(size(A,1), size(A,2),3);

for i=1:size(nedgelist, 2)
    %color = round(colors(i, :) * 255);
    Vertices=nedgelist{1,i};
    
    N0 = normal_from_gradient(Vertices);
    N  = normal_from_line(Vertices, 10);
    N  = N .* repmat((2 * (sum(N .* N0, 2) > 0) - 1), 1, 2);   % recommanded, use N0 to fix the direction


%     figure;
%     plot(nedgelist{1,i}(:,2),nedgelist{1,i}(:,1), 'r.');hold on;
%     plot([Vertices(:,2) Vertices(:,2)+10*N(:,2)]', [Vertices(:,1) Vertices(:,1)+10*N(:,1)]');
%     axis equal;
    for j=1:size(nedgelist{1,i},1)
        fprintf(fileID, '%d\t%d\t%d\t%f\t%f\n', i, nedgelist{1,i}(j,1), nedgelist{1,i}(j,2), N(j,1), N(j,2));
    end
end
fclose(fileID);

% EndPoints=zeros(size(nedgelist, 2), 2, 2);
% for i=1:size(nedgelist, 2)
%     Vertices=nedgelist{1,i};
%     EndPoints(i,1,:)=Vertices(1,:);
%     EndPoints(i,2,:)=Vertices(end,:);
% end
% 
% for i=1:size(EndPoints,1)
%     for j=1:2
%         EndPoints(i,j,:)
%     end
% end
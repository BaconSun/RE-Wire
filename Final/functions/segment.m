function nedgelist = segment(A)
    B = ones(size(A,1), size(A,2));
    %imshow(A); 
    %figure
    for i=1:size(A,1)
        for j=1:size(A,2)
            if(A(i,j)==255) 
                B(i,j)=0;
            end
        end
    end
    %im = edge(B,'canny',0.0, 0.1);

    im_ori=B;
    % im = bwmorph(im_ori,'close');


    im  = imdilate(im_ori,strel('square',8));
    % imshow(im);
    % figure
    im  = imerode(im,strel('square',8));
    %imshow(im);

    [edgelist, labelededgeim] = edgelink(im, 10);
    drawedgelist(edgelist, size(im), 1, 'rand', 5);  axis off
    seglist   = lineseg(edgelist, 2);
    drawedgelist(seglist, size(im), 2, 'rand', 6);   axis off
    minlength = 100;
    nedgelist = cleanedgelist(edgelist, minlength);
    drawedgelist(nedgelist, size(im), 1, 'rand', 7); axis off
end
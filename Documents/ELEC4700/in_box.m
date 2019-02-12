function boxnum = in_box(position, boxes)
    boxnum = 0;
    for i=1:size(boxes,1)
        if(position(1) > boxes(i,1) && position(1) < boxes(i,2) && position(2) > boxes(i,3) && position(2) < boxes(i,4))
            boxnum = i;
            return;
        end
    end
end

function [maxdistancewithpartition, pair_nodes] = maxdist(subset)

% subset contains the x-y coordinates of the points in a region. this
% algorithm split the region 3 by 3 cells and compute the maximum distance
% between pairs in an efficient way.


[coor_subset, coor_ppb] = location(subset,3);  %% 3 by 3 partition, 9 boxes!

if sum(coor_ppb(:,2)) == 0

disp("Empty Cell")
maxdistancewithpartition = 0;
pair_nodes = [];
    return
else

max_v1 = max(subset(:,1));  % max x coordinate
max_v2 = max(subset(:,2)); % max y coordinate

min_v1 = min(subset(:,1)); % min x coor
min_v2 = min(subset(:,2)); % min y coor

dist_threshold = max(max_v1 - min_v1, max_v2- min_v2); 
dist_threshold = dist_threshold/3; % grid width.

end

if dist_threshold == 0
maxdistancewithpartition = 1e-12;
pair_nodes = [];
return
end


% if coor_ppb(1,2) > 0 && coor_ppb(9,2) > 0 || coor_ppb(3,2) > 0 && coor_ppb(7,2) > 0

    dat1 = subset(coor_subset(:,4) == 1, :);
    dat9 = subset(coor_subset(:,4) == 9, :);

    dat3 = subset(coor_subset(:,4) == 3, :);
    dat7 = subset(coor_subset(:,4) == 7, :);

% md1 = dist(dat1,dat9);  %[
[md1, pts_idx1] = dist(dat1,dat9);  % dist() gives max distance between g1 and g2
% md2 = dist(dat3, dat7); 
[md2, pts_idx2] = dist(dat3, dat7);

max1 = max([md1, md2]);
% 
if max1 ~= 0 
        if max1 == md1
        pair_nodes = [dat1(pts_idx1(1),:),dat9(pts_idx1(2),:)];
        else
        pair_nodes = [dat3(pts_idx2(1),:), dat7(pts_idx2(2),:)]; 
        end
end

maxdistancewithpartition = max1;     

if max1 < sqrt(13)*dist_threshold
        
        dat2 = subset(coor_subset(:,4) == 2, :);
        dat4 = subset(coor_subset(:,4) == 4, :);
        dat5 = subset(coor_subset(:,4) == 5, :);
        dat6 = subset(coor_subset(:,4) == 6, :);
        dat7 = subset(coor_subset(:,4) == 7, :);
        dat8 = subset(coor_subset(:,4) == 8, :);
       
        [md3, pts_idx3] = dist(dat1,dat8);
        [md4, pts_idx4] = dist(dat1,dat6);
        [md5, pts_idx5] = dist(dat4,dat3);
        [md6, pts_idx6] = dist(dat4,dat9);
        [md7, pts_idx7] = dist(dat7,dat2);
        [md8, pts_idx8] = dist(dat7,dat6);
        [md9, pts_idx9] = dist(dat3,dat8);
        [md10, pts_idx10] = dist(dat2,dat9);
max2 = max([max1, md3, md4, md5, md6, md7, md8, md9, md10]);
maxdistancewithpartition = max2; 

        if max2 ~= 0
            if max2 == md3
            pair_nodes = [dat1(pts_idx3(1),:), dat8(pts_idx3(2),:)];
            elseif max2 == md4
            pair_nodes = [dat1(pts_idx4(1),:), dat6(pts_idx4(2),:)];
            elseif max2 == md5
            pair_nodes = [dat4(pts_idx5(1),:), dat3(pts_idx5(2),:)];
            elseif max2 == md6
            pair_nodes = [dat4(pts_idx6(1),:), dat9(pts_idx6(2),:)];
            elseif max2 == md7
            pair_nodes = [dat7(pts_idx7(1),:), dat2(pts_idx7(2),:)];
            elseif max2 == md8
            pair_nodes = [dat7(pts_idx8(1),:), dat6(pts_idx8(2),:)];
            elseif max2 == md9
            pair_nodes = [dat3(pts_idx9(1),:), dat8(pts_idx9(2),:)];
            elseif max2 == md10
            pair_nodes = [dat2(pts_idx10(1),:), dat9(pts_idx10(2),:)];
            end
        end
    

    if max2 < sqrt(10)*dist_threshold
               
          [md11, pts_idx11] = dist(dat3, dat9);  
          [md12, pts_idx12] = dist(dat2, dat8);
          [md13, pts_idx13] = dist(dat1, dat7);
          [md14, pts_idx14] = dist(dat1, dat3);
          [md15, pts_idx15] = dist(dat6, dat4);
          [md16, pts_idx16] = dist(dat9, dat7);
        max3 = max([max2, md11, md12, md13, md14, md15, md16]);
        maxdistancewithpartition = max3;
    
          if max3 ~= 0
            if max3 == md11
            pair_nodes = [dat3(pts_idx11(1), :), dat9(pts_idx11(2), :)];
            elseif max3 == md12
            pair_nodes = [dat2(pts_idx12(1), :), dat8(pts_idx12(2), :)];
            elseif max3 == md13
            pair_nodes = [dat1(pts_idx13(1), :), dat7(pts_idx13(2), :)];
            elseif max3 == md14
            pair_nodes = [dat1(pts_idx14(1), :), dat3(pts_idx14(2), :)];
            elseif max3 == md15
            pair_nodes = [dat6(pts_idx15(1), :), dat4(pts_idx15(2), :)];
            elseif max3 == md16
            pair_nodes = [dat9(pts_idx16(1), :), dat7(pts_idx16(2), :)];
            end
          end


                if max3 < sqrt(8)*dist_threshold
                
[md17, pts_idx17] = dist(dat3, dat5); [md18, pts_idx18] = dist(dat9, dat5); [md19, pts_idx19] = dist(dat7, dat5); [md20, pts_idx20] = dist(dat1, dat5);
[md21, pts_idx21] = dist(dat6, dat8); [md22,pts_idx22]  = dist(dat2, dat6); [md23, pts_idx23] = dist(dat2, dat4); [md24, pts_idx24] = dist(dat8, dat4);
       
                max4 = max([max3, md17, md18, md19, md20, md21, md22, md23, md24]);
                maxdistancewithpartition = max4;
                
                    if max4 ~= 0 
                        if max4 == md17
                        pair_nodes = [dat3(pts_idx17(1), :), dat5(pts_idx17(2), :)];
                        elseif max4 == md18
                        pair_nodes = [dat9(pts_idx18(1), :), dat5(pts_idx18(2), :)];
                        elseif max4 == md19
                        pair_nodes = [dat7(pts_idx19(1), :), dat5(pts_idx19(2), :)];
                        elseif max4 == md20
                        pair_nodes = [dat1(pts_idx20(1), :), dat5(pts_idx20(2), :)];
                        elseif max4 == md21
                        pair_nodes = [dat6(pts_idx21(1), :), dat8(pts_idx21(2), :)];
                        elseif max4 == md22
                        pair_nodes = [dat2(pts_idx22(1), :), dat6(pts_idx22(2), :)];
                        elseif max4 == md23
                        pair_nodes = [dat2(pts_idx23(1), :), dat4(pts_idx23(2), :)];
                        elseif max4 == md24
                        pair_nodes = [dat8(pts_idx24(1), :), dat4(pts_idx24(2), :)];
                        end
                    end


                            if max4 < sqrt(4)*dist_threshold
                         
[md25, pts_idx25] = dist(dat3, dat6); [md26, pts_idx26] = dist(dat3, dat2); [md27, pts_idx27]= dist(dat9, dat6); [md28, pts_idx28] = dist(dat5, dat6);
[md29, pts_idx29] = dist(dat9, dat8); [md30, pts_idx30] = dist(dat8, dat5); [md31,pts_idx31] = dist(dat8, dat7); [md32, pts_idx32] = dist(dat4, dat7);
[md33, pts_idx33] = dist(dat2, dat5); [md34, pts_idx34] = dist(dat5, dat4); [md35, pts_idx35]= dist(dat1, dat4); [md36, pts_idx36] = dist(dat1, dat2);
                            
        max5 = max([max4, md25, md26, md27, md28, md29, md30, md31, md32, md33, md34, md35, md36]);
        maxdistancewithpartition = max5;

        if max5~= 0
            if max5 == md25
            pair_nodes = [dat3(pts_idx25(1), :), dat6(pts_idx25(2), :)];
            elseif max5 == md26
            pair_nodes = [dat3(pts_idx26(1), :), dat2(pts_idx26(2), :)];
            elseif max5 == md27
            pair_nodes = [dat9(pts_idx27(1), :), dat6(pts_idx27(2), :)];
            elseif max5 == md28
            pair_nodes = [dat5(pts_idx28(1), :), dat6(pts_idx28(2), :)];
            elseif max5 == md29
            pair_nodes = [dat9(pts_idx29(1), :), dat8(pts_idx29(2), :)];
            elseif max5 == md30
            pair_nodes = [dat8(pts_idx30(1), :), dat5(pts_idx30(2), :)];
            elseif max5 == md31
            pair_nodes = [dat8(pts_idx31(1), :), dat7(pts_idx31(2), :)];
            elseif max5 == md32
            pair_nodes = [dat4(pts_idx32(1), :), dat7(pts_idx32(2), :)];
            elseif max5 == md33
            pair_nodes = [dat2(pts_idx33(1), :), dat5(pts_idx33(2), :)];
            elseif max5 == md34
            pair_nodes = [dat5(pts_idx34(1), :), dat4(pts_idx34(2), :)];
            elseif max5 == md35
            pair_nodes = [dat1(pts_idx35(1), :), dat4(pts_idx35(2), :)];
            elseif max5 == md36
            pair_nodes = [dat1(pts_idx36(1), :), dat2(pts_idx36(2), :)];
            end
        

             
                                     if max5 < sqrt(2)*dist_threshold
[md37, pts_idx37] = dist(dat1, dat1); [md38, pts_idx38]= dist(dat2, dat2); [md39, pts_idx39]= dist(dat3, dat3);[md40, pts_idx40]= dist(dat4, dat4); [md41, pts_idx41] = dist(dat5, dat5);
[md42, pts_idx42] = dist(dat6, dat6);[md43, pts_idx43] = dist(dat7, dat7); [md44, pts_idx44] = dist(dat8, dat8); [md45, pts_idx45] = dist(dat9, dat9);
                        max6 = max([max5, md37, md38, md39, md40, md41, md42, md43, md44, md45]);
                        maxdistancewithpartition = max6;
                                     end
        if max6~= 0
            if max6 == md37
            pair_nodes = [dat1(pts_idx37(1), :), dat1(pts_idx37(2), :)];
            elseif max6 == md38
            pair_nodes = [dat2(pts_idx38(1), :), dat2(pts_idx38(2), :)];
            elseif max6 == md39
            pair_nodes = [dat3(pts_idx39(1), :), dat3(pts_idx39(2), :)];
            elseif max6 == md40
            pair_nodes = [dat4(pts_idx40(1), :), dat4(pts_idx40(2), :)];
            elseif max6 == md41
            pair_nodes = [dat5(pts_idx41(1), :), dat5(pts_idx41(2), :)];
            elseif max6 == md42
            pair_nodes = [dat6(pts_idx42(1), :), dat6(pts_idx42(2), :)];
            elseif max6 == md43
            pair_nodes = [dat7(pts_idx43(1), :), dat7(pts_idx43(2), :)];
            elseif max6 == md44
            pair_nodes = [dat8(pts_idx44(1), :), dat8(pts_idx44(2), :)];
            elseif max6 == md45
            pair_nodes = [dat9(pts_idx45(1), :), dat9(pts_idx45(2), :)];
            end
 
        end
        end                    
        end
                end
        end
end
end



function [max_distance,pts_idx] = dist(dat1,dat2)

l1 = size(dat1,1);
l2 = size(dat2,1);

if l1 ==0 || l2 ==0
    max_distance = 0;
    pts_idx = [];
    return
else
    distance = zeros(l1*l2,3);
    counter = 0;
        for i=1:l1
            for j=1:l2
                counter = counter + 1;
                 distance(counter,1)=i;
                 distance(counter,2)=j;
                distance(counter,3) = sqrt((dat1(i,1)-dat2(j,1))^2 + (dat1(i,2)-dat2(j,2))^2);
            end
        end

distance = sortrows(distance, 3, 'descend');
max_distance = distance(1,3);
pts_idx = distance(1,1:2);
end


end


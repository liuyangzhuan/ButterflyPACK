clear all;
clc;
close all;
edgecen = load('fort.113');
Maxlevel = max(edgecen(:,1))
for level=0:Maxlevel
    figure(level+1)
    edgecen_onelevel = edgecen(find(edgecen(:,1)==level),:);
    for group = 2^level:2^(level+1)-1
        edgecen_onegroup = edgecen_onelevel(find(edgecen_onelevel(:,2)==group),:);
        plot3(edgecen_onegroup(:,3),edgecen_onegroup(:,4),edgecen_onegroup(:,5),'.','Color', [rand, rand, rand])
        hold on        
    end
    axis equal
end

% skel = load('fort.114');
% figure(2)
% plot3(skel(:,1),skel(:,2),skel(:,3),'o','Color', 'k')
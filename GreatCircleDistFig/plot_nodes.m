% For a 4 deg RLL mesh ix = 3750
function plot_nodes(nodefile, facefile, ix, dist, disttype)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes
nodes = load(nodefile);

%plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.');
hold on;
%nn = size(nodes,1);
%plot3(nodes(83:nn,1), nodes(83:nn,2), nodes(83:nn,3), 'r.');
gix=[]+1;
for i=1:length(gix)
    plot3(nodes(gix(i),1), nodes(gix(i),2), nodes(gix(i),3), 'bo');
end
hold off;

grid;

% Plot faces
faces = load(facefile);

xnodes = nodes(:,1);
ynodes = nodes(:,2);
znodes = nodes(:,3);

xfaces = xnodes(faces(:,:));
yfaces = ynodes(faces(:,:));
zfaces = znodes(faces(:,:));

if(size(faces,1) == 1)
    xfaces = xfaces.';
    yfaces = yfaces.';
    zfaces = zfaces.';
end

hold on;
axis off;

x0 = sum(xfaces(ix,:))/length(xfaces(ix,:));
y0 = sum(yfaces(ix,:))/length(yfaces(ix,:));
z0 = sum(zfaces(ix,:))/length(zfaces(ix,:));

lon0 = atan2(y0, x0);
lat0 = asin(z0 / sqrt(x0^2 + y0^2 + z0^2));

[lon0*180/pi lat0*180/pi]

for i = 1:size(xfaces,1)
    if (ix == i)
        patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),[0.3 0.3 0.3],'EdgeColor','k','LineWidth',2);
    else
        x1 = sum(xfaces(i,:))/length(xfaces(i,:));
        y1 = sum(yfaces(i,:))/length(yfaces(i,:));
        z1 = sum(zfaces(i,:))/length(zfaces(i,:));
        
        lon1 = atan2(y1, x1);
        lat1 = asin(z1 / sqrt(x1^2 + y1^2 + z1^2));
        
        gcd = acos(sin(lat1)*sin(lat0) + cos(lat1)*cos(lat0)*cos(lon1-lon0));
        lld = sqrt((lat1-lat0)^2 + (lon1-lon0)^2);
        
        if ((disttype == 0) && (lld < dist))
            patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),[0.6 0.6 0.6],'EdgeColor','k','LineWidth',2);
        elseif ((disttype == 1) && (gcd < dist))
            patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),[0.6 0.6 0.6],'EdgeColor','k','LineWidth',2);
        else
            patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w','EdgeColor','k','LineWidth',2);
        end
    end
end
%patch(xfaces(1,:), yfaces(1,:), zfaces(1,:),'g','EdgeColor','k','LineWidth',2);
%patch(xfaces(690,:), yfaces(690,:), zfaces(690,:),'b','EdgeColor','k','LineWidth',2);
%patch(xfaces(1336,:), yfaces(1336,:), zfaces(1336,:),'r','EdgeColor','k','LineWidth',2);
%patch(xfaces(11816,:), yfaces(11816,:), zfaces(11816,:),'r');
hold off;

view([-45 60]);
camzoom(3.0);

set(gca, 'FontSize', 14);
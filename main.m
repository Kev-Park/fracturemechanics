%% Using n-CST mesh
clear;

nodes = load("side_circle_nodes.dat");
nodes = nodes(:,2:3);

elements = load("side_circle_elements.dat");
elements = elements(:,6:8);

E = 10^11;
nu = .1;
h = 0.1;

% Constitutive matrix
C = (E/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% Calculate area, B-matrix, K-matrix (element-wise)
for i = 1:length(elements)
    x1 = nodes(elements(i,1),1);
    y1 = nodes(elements(i,1),2);
    x2 = nodes(elements(i,2),1);
    y2 = nodes(elements(i,2),2);
    x3 = nodes(elements(i,3),1);
    y3 = nodes(elements(i,3),2);

    A_{i} = 0.5 * det([1 x1 y1; 1 x2 y2; 1 x3 y3]);

    B_{i} = (1/(2*A_{i})) * [y2-y3 0 y3-y1 0 y1-y2 0; 0 x3-x2 0 x1-x3 0 x2-x1; x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];

    Ke_{i} = h*A_{i}*B_{i}.'*C*B_{i}; 
end

K_global = zeros(2*length(nodes),2*length(nodes));

% Assemble global K-matrix
for i = 1:length(elements)
    
    ids = [2*elements(i,1) - 1; 2*elements(i,1); 2*elements(i,2) - 1; 2*elements(i,2); 2*elements(i,3) - 1; 2*elements(i,3)];

    

    for j = 1:length(Ke_{i}(:,1)) % Rows
       
        for k = 1:length(Ke_{i}(1,:)) % Columns



            K_global(ids(j),ids(k)) = Ke_{i}(j,k) + K_global(ids(j),ids(k));

        end

    end

end

% Assemble F_ext and enforce boundary conditions
% F_ext = zeros(2*length(nodes),1);
% 
% F = 0.5;
% K_global_enforced = K_global;
% 
% for i = 1:length(nodes)
% 
%     if nodes(i,2) == 10 % If y is at top, apply stress
%         F_ext(2*i) = -F; % Apply to y-dir
%     end
% 
%     if nodes(i,2) == 0 % If y is at bottom, apply stress
%         F_ext(2*i) = F; % Apply to y-dir
%     end
% 
%     if nodes(i,1) == 10 % If x-coordinate is along the right boundary (fix)
%         K_global_enforced(2*i-1,:) = 0;
%         K_global_enforced(2*i-1,2*i-1) = 1;
%         K_global_enforced(2*i,:) = 0;
%         K_global_enforced(2*i,2*i) = 1;
%     end
% end

num_nodes = length(nodes);
nelements = length(elements);
dim = num_nodes * 2;
F = 30000;
H = 10;
W = 10;

top_right_node_num = 0;
bot_right_node_num = 0;
top_dof_forced = [];
bot_dof_forced = [];
for i = 1:num_nodes
    if (nodes(i, 2) == H && nodes(i, 1) == W)
        top_right_node_num = i;
    elseif (nodes(i, 2) == 0 && nodes(i, 1) == W)
        bot_right_node_num = i;
    end

    if (nodes(i, 2) == H)
        top_dof_forced = cat(2, top_dof_forced, [i*2-1, i*2]);
    elseif (nodes(i, 2) == 0)
        bot_dof_forced = cat(2, bot_dof_forced, [i*2-1, i*2]);
    end
end

dof_reactive = [top_right_node_num*2-1, bot_right_node_num*2-1, bot_right_node_num*2];

dof_active = zeros(dim-length(dof_reactive), 1);

activeIndex = 1;
for i = 1:dim
    index = find(dof_reactive==i);
    if isempty(index)
        dof_active(activeIndex) = i;
        activeIndex = activeIndex + 1;
    end
end


top_F_per_node = F/(length(top_dof_forced)/2);
bot_F_per_node = -F/(length(bot_dof_forced)/2);

F = zeros(dim, 1);
for i = 1:dim
    index = find(top_dof_forced==i);
    if ~isempty(index) && (-1)^index == 1
        F(i, 1) = top_F_per_node;
    end

    index2 = find(bot_dof_forced==i);
    if ~isempty(index2) && (-1)^index2 == 1
        F(i, 1) = bot_F_per_node;
    end
end

% Enforcing Boundary Conditions in K_global, F

BCs = zeros(dim, 1); % All zeros due to no set displacements

K_global_enforced = K_global;
F_ext = F;

K_global_enforced(top_right_node_num*2-1, :) = zeros(1, dim);
K_global_enforced(top_right_node_num*2-1, top_right_node_num*2-1) = 1;
F_ext(top_right_node_num*2-1, 1) = BCs(i, 1);

K_global_enforced(top_right_node_num*2, :) = zeros(1, dim);
K_global_enforced(top_right_node_num*2, top_right_node_num*2) = 1;
F_ext(top_right_node_num*2, 1) = BCs(i, 1);

K_global_enforced(bot_right_node_num*2-1, :) = zeros(1, dim);
K_global_enforced(bot_right_node_num*2-1, bot_right_node_num*2-1) = 1;
F_ext(bot_right_node_num*2-1, 1) = BCs(bot_right_node_num*2-1, 1);

K_global_enforced(bot_right_node_num*2, :) = zeros(1, dim);
K_global_enforced(bot_right_node_num*2, bot_right_node_num*2) = 1;
F_ext(bot_right_node_num*2, 1) = BCs(bot_right_node_num*2, 1);



% Find displacement
d = K_global_enforced\F_ext;

% Find stresses
% sigma_vm = zeros(length(elements),1);
% for i = 1:length(elements)
%     node1 = elements(i,1);
%     node2 = elements(i,2);
%     node3 = elements(i,3);
%     d_i = [nodes(node1, 1) nodes(node1, 2) nodes(node2, 1) nodes(node2, 2) nodes(node3, 1) nodes(node3, 2)]';
% 
%     sigma = C*B_{i}*d_i;
% 
%     sigma_vm(i) = sqrt(sigma(1)^2 - sigma(1)*sigma(2) + sigma(2)^2 + 3*sigma(3)^2);
% end

sigma_vm = zeros(nelements ,1);
sigmas = zeros(nelements,3);
for i = 1:nelements
    node1 = elements(i,1);
    node2 = elements(i,2);
    node3 = elements(i,3);
    % node_pos = cat(1, nodes(node1, :), nodes(node2, :), nodes(node3, :));

    d1 = [d(node1*2-1), d(node1*2)];
    d2 = [d(node2*2-1), d(node2*2)];
    d3 = [d(node3*2-1), d(node3*2)];
    d_i = cat(2, d1, d2, d3);
    d_i = reshape(d_i, [6, 1]);

    sigma = C*B_{i}*d_i;

    sigmas(i,:) = sigma;
    sigma_vm(i) = sqrt(sigma(1)^2 - sigma(1)*sigma(2) + sigma(2)^2 + (3/4)*sigma(3)^2);
end

% Displaced nodes
nodes_displaced = nodes + reshape(d, 2, [])';

% Plot triangles

% Apply scaling to von Mises stresses so coloration is more apparent
sigma_vm_log = log10(sigma_vm);

patch('Faces',elements,'Vertices',nodes_displaced,'FaceVertexCData', sigma_vm_log,'FaceColor','flat','EdgeColor','none', 'ButtonDownFcn', @(src, event) showStress(src, event, sigmas));
colormap('jet');

% Find tick marks
ticks_real = logspace(log10(min(sigma_vm)), log10(max(sigma_vm)), 10);
ticks_log = log10(ticks_real);

colorbar('Ticks', ticks_log, 'TickLabels', ticks_real);

% Display stress on hover
function showStress(src, event, sigmas)

    pt = event.IntersectionPoint(1:2);
    faces = src.Faces;
    verts = src.Vertices;

    for i = 1:size(faces,1)
        tri_idx = faces(i,:);
        tri = verts(tri_idx, 1:2);  % 2D vertices

        if inpolygon(pt(1), pt(2), tri(:,1), tri(:,2))
            % fea
            sigmas(i,:)

            % analytical - hard coded values for formula (!)
            K = 30000 * sqrt(pi * 3) * (1.122 - 0.231 * (3/10) + 10.55 * (3/10)^2 - 21.71 * (3/10)^3 + 30.382 * (3/10)^4);
            r = sqrt((pt(1)-3)^2 + (pt(2)-5)^2);
            theta = atan((pt(2)-5)/(pt(1)-3));
            xx = (K/sqrt(2*pi*r)) * cos(theta/2) * (1 - sin(theta/2) * sin(3*theta/2));
            yy = (K/sqrt(2*pi*r)) * cos(theta/2) * (1 + sin(theta/2) * sin(3*theta/2));
            xy = (K/sqrt(2*pi*r)) * cos(theta/2) * sin(theta/2) * cos(3*theta/2);

            % uncomment out as desired
            disp("analytic");
            [xx, yy, 2*xy]
            r
            theta

            return;
        end
    end
end


% plot(nodes(:,1),nodes(:,2),'.');
% max(d)
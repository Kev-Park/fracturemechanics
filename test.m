%% Function to assemble stiffness matrix, stress, and strain

% Assuming plane stress
function [Ke, strain, stress] = calcCST(node_pos, d, E, v, h)
    x1 = node_pos(1, 1);
    x2 = node_pos(2, 1);
    x3 = node_pos(3, 1);
    y1 = node_pos(1, 2);
    y2 = node_pos(2, 2);
    y3 = node_pos(3, 2);
    
    J = 0.5*[x2-x1 y2-y1; x3-x1 y3-y1];
    
    A1 = [1 0 0 0; 0 0 0 1; 0 1 1 0];
    A2 = cat(2, cat(1, inv(J), zeros(2, 2)), cat(1, zeros(2, 2), inv(J)));
    A3 = 0.5*[-1 0 1 0 0 0; -1 0 0 0 1 0; 0 -1 0 1 0 0; 0 -1 0 0 0 1];
    
    B = A1*A2*A3;
    
    C = E / (1-v^2) * [1 v 0; v 1 0; 0 0 (1-v)/2];
    
    strain = B*d;
    stress = C*strain;
    Ke = h*2*det(J)*B.'*C*B;

end

%% Function to assemble global stiffness matrix

function Kg = assembleKg(K_element_all, num_nodes)
    nelements = size(K_element_all, 1);
    dof_per_element = size(K_element_all{1, 1}, 1);
    dim = num_nodes * 2;
    Kg = zeros(dim, dim);

    for i = 1:nelements
    
        I = zeros(3, 1);
        I(1) = K_element_all{i, 2}(1) * 2 - 1;
        I(2) = K_element_all{i, 2}(1) * 2;
        I(3) = K_element_all{i, 2}(2) * 2 - 1;
        I(4) = K_element_all{i, 2}(2) * 2;
        I(5) = K_element_all{i, 2}(3) * 2 - 1;
        I(6) = K_element_all{i, 2}(3) * 2;
    
        for rowIndex = 1:dof_per_element
            for colIndex = 1:dof_per_element
                ii = I(rowIndex);
                jj = I(colIndex);
                Kg(ii, jj) = Kg(ii, jj) + K_element_all{i, 1}(rowIndex, colIndex);
            end
        end
    
    end
end


%% Running

% For more nodes
elements = load('elements.dat');
nodes = load('nodes.dat');
nelements = length(elements);
num_nodes = length(nodes);
node_positions = nodes(:, 2:3);
e_node_nums = elements(:, 6:8);


K_element_all = cell(nelements, 2);

d = zeros(6, 1); % Just set to 0 for an input, not used
E = 210000000000; % Pa
v = 0.2;
h = 0.1; % m
F = 3000; % N
W = 10;
H = 10;

for i = 1:nelements
    node_pos = cat(1, node_positions(e_node_nums(i, 1), :), node_positions(e_node_nums(i, 2), :), node_positions(e_node_nums(i, 3), :));
    [Ke, ~, ~] = calcCST(node_pos, d, E, v, h);

    K_element_all{i, 1} = Ke;
    K_element_all{i, 2} = e_node_nums(i, :);
end

Kg = assembleKg(K_element_all, num_nodes);



% Calculating the element nodal force vector

dim = num_nodes * 2;

top_right_node_num = 0;
bot_right_node_num = 0;
top_dof_forced = [];
bot_dof_forced = [];
for i = 1:num_nodes
    if (node_positions(i, 2) == H && node_positions(i, 1) == W)
        top_right_node_num = i;
    elseif (node_positions(i, 2) == 0 && node_positions(i, 1) == W)
        bot_right_node_num = i;
    end

    if (node_positions(i, 2) == H)
        top_dof_forced = cat(2, top_dof_forced, [i*2-1, i*2]);
    elseif (node_positions(i, 2) == 0)
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

% Enforcing Boundary Conditions in Kg, F

BCs = zeros(dim, 1); % All zeros due to no set displacements

Kg_BC = Kg;
F_BC = F;

Kg_BC(top_right_node_num*2-1, :) = zeros(1, dim);
Kg_BC(top_right_node_num*2-1, top_right_node_num*2-1) = 1;
F_BC(top_right_node_num*2-1, 1) = BCs(i, 1);

Kg_BC(top_right_node_num*2, :) = zeros(1, dim);
Kg_BC(top_right_node_num*2, top_right_node_num*2) = 1;
F_BC(top_right_node_num*2, 1) = BCs(i, 1);

Kg_BC(bot_right_node_num*2-1, :) = zeros(1, dim);
Kg_BC(bot_right_node_num*2-1, bot_right_node_num*2-1) = 1;
F_BC(bot_right_node_num*2-1, 1) = BCs(bot_right_node_num*2-1, 1);

Kg_BC(bot_right_node_num*2, :) = zeros(1, dim);
Kg_BC(bot_right_node_num*2, bot_right_node_num*2) = 1;
F_BC(bot_right_node_num*2, 1) = BCs(bot_right_node_num*2, 1);


% Solving for displacements

du = Kg_BC\F_BC; % in m

% Find stresses
sigma_vm = zeros(nelements ,1);
for i = 1:nelements
    node1 = e_node_nums(i,1);
    node2 = e_node_nums(i,2);
    node3 = e_node_nums(i,3);
    node_pos = cat(1, node_positions(node1, :), node_positions(node2, :), node_positions(node3, :));

    d1 = [du(node1*2-1), du(node1*2)];
    d2 = [du(node2*2-1), du(node2*2)];
    d3 = [du(node3*2-1), du(node3*2)];
    d = cat(2, d1, d2, d3);
    d = reshape(d, [6, 1]);

    [~, sigma, ~] = calcCST(node_pos, d, E, v, h);

    sigma_vm(i) = sqrt(sigma(1)^2 - sigma(1)*sigma(2) + sigma(2)^2 + 3*sigma(3)^2);
end

%% Plot triangles

% Displaced nodes
nodes_displaced = node_positions + reshape(du, 2, [])';

% Plot triangles

% Apply scaling to von Mises stresses so coloration is more apparent
sigma_vm_sq = sigma_vm.^1;

patch('Faces',e_node_nums,'Vertices',nodes_displaced,'FaceVertexCData', sigma_vm_sq,'FaceColor','flat','EdgeColor','none');
colormap('hot');
colorbar;
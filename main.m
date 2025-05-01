%% Using n-CST mesh
clear;

nodes = load("nodes.dat");
nodes = nodes(:,2:3);

elements = load("elements.dat");
elements = elements(:,6:8);

E = 10^11;
nu = 0.0;
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
F_ext = zeros(2*length(nodes),1);

F = .1;
K_global_enforced = K_global;

for i = 1:length(nodes)

    if nodes(i,2) == 10 % If y is at top, apply stress
        F_ext(2*i) = -F; % Apply to y-dir
    end

    if nodes(i,2) == 0 % If y is at bottom, apply stress
        F_ext(2*i) = F; % Apply to y-dir
    end

    if nodes(i,1) == 10 % If x-coordinate is along the right boundary (fix)
        K_global_enforced(2*i-1,:) = 0;
        K_global_enforced(2*i-1,2*i-1) = 1;
        K_global_enforced(2*i,:) = 0;
        K_global_enforced(2*i,2*i) = 1;
    end
end

% Find displacement
d = inv(K_global_enforced)*F_ext;

% Displaced nodes
nodes_displaced = nodes + reshape(d, 2, [])';

% Plot triangles
patch('Faces',elements,'Vertices',nodes_displaced,'FaceColor','cyan','EdgeColor','black');

% plot(nodes(:,1),nodes(:,2),'.');
% max(d)
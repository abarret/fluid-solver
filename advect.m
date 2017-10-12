function vals = advect(q,u,v)
% Computes u*grad(q) using a wave propagation method
k = 3; gcw = k;
global dx; global dy;
current_vals = fillBoundariesCenter(q, gcw);
% Reconstruct values
[left_vals_x, right_vals_x, low_vals_y, up_vals_y] = WENO_reconstruction_2D(current_vals, gcw, k);
% Calculate waves in x
wave_left_x = max(u(:,1:end-1),0.0).*(left_vals_x(:,1:end-1)-right_vals_x(:,1:end-1));
wave_right_x = min(u(:,2:end),0.0).*(left_vals_x(:,2:end)-right_vals_x(:,2:end));
wave_center_x = (u(:,1:end-1)+u(:,2:end))*0.5.*(right_vals_x(:,2:end) - left_vals_x(:,1:end-1));
% Calculate waves in y
wave_left_y = max(v(1:end-1,:),0.0).*(low_vals_y(1:end-1,:)-up_vals_y(1:end-1,:));
wave_right_y = min(v(2:end,:),0.0).*(low_vals_y(2:end,:)-up_vals_y(2:end,:));
wave_center_y = (v(1:end-1,:)+v(2:end,:))*0.5.*(up_vals_y(2:end,:) - low_vals_y(1:end-1,:));

vals = -1.0/dx*(wave_left_x + wave_right_x + wave_center_x) - 1.0/dy*(wave_left_y + wave_right_y + wave_center_y);
vals = vals;
end

function [left_vals_x, right_vals_x, low_vals_y, up_vals_y] = ...
    WENO_reconstruction_2D(cell_average, gcw, k)

global interp_weights_at_cell_sides;
global smooth_constants;
interp_weights_at_cell_sides = CalculateInterpolationCoefficients(k);
smooth_constants = CalculateConstantsForSmoothAveraging(k);

interp_values_x_right = cell(1,k);
interp_values_x_left = cell(1,k);

interp_values_y_up = cell(1,k);
interp_values_y_low = cell(1,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form interpolants at cell sides %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 0:k-1
%     interp_values_x_right{r+1} = zeros(size(cell_average(gcw:end-gcw,gcw+1:end-gcw)));
%     interp_values_x_left{r+1}  = zeros(size(cell_average(gcw+1:end-gcw+1,gcw+1:end-gcw)));
%     
%     interp_values_y_up{r+1} = zeros(size(cell_average(gcw+1:end-gcw,gcw:end-gcw)));
%     interp_values_y_low{r+1}  = zeros(size(cell_average(gcw+1:end-gcw,gcw+1:end-gcw+1)));
%     for j = 0:k-1
%         interp_values_x_right{r+1} = interp_values_x_right{r+1} + interp_weights_at_cell_sides(r+2,j+1)*cell_average(gcw-r+j:end-gcw-r+j,gcw+1:end-gcw);
%         interp_values_x_left{r+1}  = interp_values_x_left{r+1} + interp_weights_at_cell_sides(r+1,j+1)*cell_average(gcw+1-r+j:end-gcw+1-r+j,gcw+1:end-gcw);
%         
%         interp_values_y_up{r+1}   = interp_values_y_up{r+1} + interp_weights_at_cell_sides(r+2,j+1)*cell_average(gcw+1:end-gcw,gcw-r+j:end-gcw-r+j);
%         interp_values_y_low{r+1}  = interp_values_y_low{r+1} + interp_weights_at_cell_sides(r+1,j+1)*cell_average(gcw+1:end-gcw,gcw+1-r+j:end-gcw+1-r+j);
%     end
    interp_values_x_right{r+1} = zeros(size(cell_average(gcw+1:end-gcw,gcw:end-gcw)));
    interp_values_x_left{r+1}  = zeros(size(cell_average(gcw+1:end-gcw,gcw+1:end-gcw+1)));
    
    interp_values_y_up{r+1} = zeros(size(cell_average(gcw:end-gcw,gcw+1:end-gcw)));
    interp_values_y_low{r+1}  = zeros(size(cell_average(gcw+1:end-gcw+1,gcw+1:end-gcw)));
    for j = 0:k-1
        interp_values_x_right{r+1} = interp_values_x_right{r+1} + interp_weights_at_cell_sides(r+2,j+1)*cell_average(gcw+1:end-gcw,gcw-r+j:end-gcw-r+j);
        interp_values_x_left{r+1}  = interp_values_x_left{r+1} + interp_weights_at_cell_sides(r+1,j+1)*cell_average(gcw+1:end-gcw,gcw+1-r+j:end-gcw+1-r+j);
        
        interp_values_y_up{r+1}   = interp_values_y_up{r+1} + interp_weights_at_cell_sides(r+2,j+1)*cell_average(gcw-r+j:end-gcw-r+j,gcw+1:end-gcw);
        interp_values_y_low{r+1}  = interp_values_y_low{r+1} + interp_weights_at_cell_sides(r+1,j+1)*cell_average(gcw+1-r+j:end-gcw+1-r+j,gcw+1:end-gcw);
    end
end

%%%%%%%%%%%%%%%%%%%%%
% Calculate Weights %
%%%%%%%%%%%%%%%%%%%%%


smooth_id_x_right = CalculateSmoothnessIndicators(k, cell_average(gcw+1:end-gcw,gcw-(k-1):end-gcw+(k-1)), 2, k-1);
weights_x_right = CalculateWeights(k, smooth_id_x_right, smooth_constants);
smooth_id_x_left = CalculateSmoothnessIndicators(k, cell_average(gcw+1:end-gcw,gcw+1-(k-1):end-gcw+1+(k-1)), 2, k-1);
weights_x_left = CalculateWeights(k, smooth_id_x_left, smooth_constants(end:-1:1));

smooth_id_y_up = CalculateSmoothnessIndicators(k, cell_average(gcw-(k-1):end-gcw+(k-1),gcw+1:end-gcw), 1, k-1);
weights_y_up = CalculateWeights(k, smooth_id_y_up, smooth_constants);
smooth_id_y_low = CalculateSmoothnessIndicators(k, cell_average(gcw+1-(k-1):end-gcw+1+(k-1),gcw+1:end-gcw), 1, k-1);
weights_y_low = CalculateWeights(k, smooth_id_y_low, smooth_constants(end:-1:1));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Interpolant %
%%%%%%%%%%%%%%%%%%%%%%%%%
left_vals_x = zeros(size(cell_average(gcw+1:end-gcw,gcw+1:end-gcw+1)));
right_vals_x = zeros(size(cell_average(gcw+1:end-gcw,gcw:end-gcw)));
low_vals_y = zeros(size(cell_average(gcw+1:end-gcw+1,gcw+1:end-gcw)));
up_vals_y = zeros(size(cell_average(gcw:end-gcw,gcw+1:end-gcw)));
for r = 0:k-1
    right_vals_x = right_vals_x + weights_x_right{r+1}.*interp_values_x_right{r+1};
    left_vals_x = left_vals_x + weights_x_left{r+1}.*interp_values_x_left{r+1};
    
    up_vals_y = up_vals_y + weights_y_up{r+1}.*interp_values_y_up{r+1};
    low_vals_y = low_vals_y + weights_y_low{r+1}.*interp_values_y_low{r+1};
end

function weights = CalculateWeights(k, smooth_id, exact_weights)
    eps = 1.0e-6;
    weights = cell(1,k);
    tot = zeros(size(smooth_id{1}));
    for i = 1:k
        alpha = exact_weights(i)./(eps+smooth_id{i}).^2;
        weights{i} = alpha;
        tot = tot + alpha;
    end
    for i = 1:k
        weights{i} = weights{i}./tot;
    end
end
function smooth_id = CalculateSmoothnessIndicators(k, vals, dim, gcw)
    smooth_id = cell(1,k);
    if(dim == 1)
        % x direction
        switch(k)
            case 3
                smooth_id{1} = 13.0/12.0*(vals(gcw+1:end-gcw,:) - 2*vals(gcw+2:end-gcw+1,:) + vals(gcw+3:end-gcw+2,:)).^2 + ...
                    0.25*(3*vals(gcw+1:end-gcw,:) - 4*vals(gcw+2:end-gcw+1,:) + vals(gcw+3:end-gcw+2,:)).^2;
                smooth_id{2} = 13.0/12.0*(vals(gcw:end-gcw-1,:) - 2*vals(gcw+1:end-gcw,:) + vals(gcw+2:end-gcw+1,:)).^2 + ...
                    0.25*(vals(gcw:end-gcw-1,:) - vals(gcw+2:end-gcw+1,:)).^2;
                smooth_id{3} = 13.0/12.0*(vals(gcw-1:end-gcw-2,:) - 2*vals(gcw:end-gcw-1,:) + vals(gcw+1:end-gcw,:)).^2 + ...
                    0.25*(vals(gcw-1:end-gcw-2,:) - 4*vals(gcw:end-gcw-1,:) + 3*vals(gcw+1:end-gcw,:)).^2;
        end
    elseif(dim == 2)
        % y direction
        switch(k)
            case 3
                smooth_id{1} = 13.0/12.0*(vals(:,gcw+1:end-gcw) - 2*vals(:,gcw+2:end-gcw+1) + vals(:,gcw+3:end-gcw+2)).^2 + ...
                    0.25*(3*vals(:,gcw+1:end-gcw) - 4*vals(:,gcw+2:end-gcw+1) + vals(:,gcw+3:end-gcw+2)).^2;
                smooth_id{2} = 13.0/12.0*(vals(:,gcw:end-gcw-1) - 2*vals(:,gcw+1:end-gcw) + vals(:,gcw+2:end-gcw+1)).^2 + ...
                    0.25*(vals(:,gcw:end-gcw-1) - vals(:,gcw+2:end-gcw+1)).^2;
                smooth_id{3} = 13.0/12.0*(vals(:,gcw-1:end-gcw-2) - 2*vals(:,gcw:end-gcw-1) + vals(:,gcw+1:end-gcw)).^2 + ...
                    0.25*(vals(:,gcw-1:end-gcw-2) - 4*vals(:,gcw:end-gcw-1) + 3*vals(:,gcw+1:end-gcw)).^2;
        end
    end
end
end
function weights = CalculateInterpolationCoefficients(k)
    % weights is an k+1 by k matrix
    % r runs in rows, 
    % interpolant coefficients runs in columns
    weights = zeros(k+1,k);
    for r = -1:k-1
        for j = 0:k-1
            weights(r+2,j+1) = 0.0;
            prod = 1.0;
            for l=0:k-1
                if l ~= j
                prod = prod*(r-l+0.5)/(j-l);
                end
            end
            weights(r+2,j+1) = prod;
        end
    end
end
function weights = CalculateConstantsForSmoothAveraging(k)
    weights = zeros(k,1);
    switch(k)
        case 1
            weights(1) = 1.0;
        case 2
            weights(1) = 2/3;
            weights(2) = 1/3;
        case 3
            weights(1) = 5.0/16.0;
            weights(2) = 5.0/8.0;
            weights(3) = 1.0/16.0;
        otherwise
            error('Value of k not supported');
    end
end

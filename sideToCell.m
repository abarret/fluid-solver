function cell_data = sideToCell(sideX,sideY)
% Converts cell data to side data using linear interpolation
[~,Ny] = size(sideY);
[Nx,~] = size(sideX);
cell_data = zeros(Nx,Ny,2);
cell_data(:,:,1) = 0.5*(sideX(:,2:end)+sideX(:,1:end-1));
cell_data(:,:,2) = 0.5*(sideY(2:end,:)+sideY(1:end-1,:));
end

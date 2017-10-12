function [sideX, sideY] = cellToSide(cell_data)
% Converts cell data to side data using linear interpolation
gcw = 1;
cell_data = fillBoundariesCenter(cell_data,gcw);
sideY = 0.5*(cell_data(gcw:end-gcw,gcw+1:end-gcw)+cell_data(gcw+1:end-gcw+1,gcw+1:end-gcw));
sideX = 0.5*(cell_data(gcw+1:end-gcw,gcw:end-gcw)+cell_data(gcw+1:end-gcw,gcw+1:end-gcw+1));
end

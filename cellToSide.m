function [sideX, sideY] = cellToSide(cell_data)
% Converts cell data to side data using linear interpolation. Note that
% this assumes cell_data is scalar valued and interpolates the scalar value
% to each side.
% Inputs:
%   cell_data : Cell centered data to be interpolated
% Outputs:
%   sideX : Side centered x data
%   sideY : Side centered y data
gcw = 1;
cell_data = fillBoundariesCenter(cell_data,gcw);
sideY = 0.5*(cell_data(gcw:end-gcw,gcw+1:end-gcw)+cell_data(gcw+1:end-gcw+1,gcw+1:end-gcw));
sideX = 0.5*(cell_data(gcw+1:end-gcw,gcw:end-gcw)+cell_data(gcw+1:end-gcw,gcw+1:end-gcw+1));
end

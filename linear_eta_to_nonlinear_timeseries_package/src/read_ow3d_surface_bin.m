function [X, Y, eta_field, phi_field] = read_ow3d_surface_bin(filename)
%READ_OW3D_SURFACE_BIN Read an OW3D EP_XXXXX.bin surface file.

byteorder = 'ieee-le';
fid = fopen(filename, 'r', byteorder);
if fid < 0
    error('Unable to open OW3D bin file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

fread(fid, 1, 'int32');
Nx = fread(fid, 1, 'int32');
Ny = fread(fid, 1, 'int32');
fread(fid, 1, 'int32');

fread(fid, 1, 'int32');
X = fread(fid, [Nx Ny], 'float64');
Y = fread(fid, [Nx Ny], 'float64');
fread(fid, 1, 'int32');

fread(fid, 1, 'int32');
eta_field = fread(fid, [Nx Ny], 'float64');
phi_field = fread(fid, [Nx Ny], 'float64');
end

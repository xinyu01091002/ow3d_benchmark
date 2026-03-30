function data = read_ow3d_kinematics_file(file_path, phit_mode)
%READ_OW3D_KINEMATICS_FILE Read an OW3D KinematicsXX.bin file.

if nargin < 2 || isempty(phit_mode)
    phit_mode = 'uncorrected';
end

[it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, h, sigma, t] = ...
    read_kinematics_file_local_package(file_path, phit_mode); %#ok<ASGLU>

data = struct();
data.it = it;
data.eta = eta;
data.etat_m = etat_m;
data.etatt_m = etatt_m;
data.phi = phi;
data.phit = phit_m;
data.p = p_m;
data.ut = ut_m;
data.u = u;
data.v = v;
data.w = w;
data.uz = uz;
data.vz = vz;
data.wz = wz;
data.x = x;
data.y = y;
data.h = h;
data.sigma = sigma;
data.t = t;
end

function [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, h, sigma, t] = read_kinematics_file_local_package(file_path, phit_mode)
    nbits = 32;
    compute_derivatives = true;

    if nbits == 32
        int_nbit = 'int';
    elseif nbits == 64
        int_nbit = 'int64';
    else
        error('Illegal value for nbits: %d', nbits);
    end

    fid = fopen(file_path, 'r', 'ieee-le');
    if fid < 0
        error('Could not open kinematics file: %s', file_path);
    end
    cleanup = onCleanup(@() fclose(fid));

    fread(fid, 1, int_nbit);
    xbeg = fread(fid, 1, 'int');
    xend = fread(fid, 1, 'int');
    xstride = fread(fid, 1, 'int');
    ybeg = fread(fid, 1, 'int');
    yend = fread(fid, 1, 'int');
    ystride = fread(fid, 1, 'int');
    tbeg = fread(fid, 1, 'int');
    tend = fread(fid, 1, 'int');
    tstride = fread(fid, 1, 'int');
    dt = fread(fid, 1, 'double');
    nz = fread(fid, 1, 'int');
    sigma = zeros(nz, 1);
    fread(fid, 2, int_nbit);

    nx = floor((xend - xbeg) / xstride) + 1;
    ny = floor((yend - ybeg) / ystride) + 1;
    nt = floor((tend - tbeg) / tstride) + 1;

    tmp = zeros(nx * ny * max(nz, 5), 1);
    tmp(1:5 * nx * ny) = fread(fid, 5 * nx * ny, 'double');
    fread(fid, 2, int_nbit);

    x = zeros(nx, ny);
    y = zeros(nx, ny);
    h = zeros(nx, ny);
    x(:) = tmp(1:5:5 * nx * ny);
    y(:) = tmp(2:5:5 * nx * ny);
    h(:) = tmp(3:5:5 * nx * ny);

    for i = 1:nz
        sigma(i) = fread(fid, 1, 'double');
    end
    fread(fid, 2, int_nbit);

    eta = zeros(nt, nx, ny);
    etax = zeros(nt, nx, ny);
    etay = zeros(nt, nx, ny);
    phi = zeros(nt, nz, nx, ny);
    w = zeros(nt, nz, nx, ny);
    u = zeros(nt, nz, nx, ny);
    uz = zeros(nt, nz, nx, ny);
    v = zeros(nt, nz, nx, ny);
    vz = zeros(nt, nz, nx, ny);
    wz = zeros(nt, nz, nx, ny);
    t = (0:nt - 1) * dt * tstride;

    it = 0;
    for it_read = 1:nt - 1
        tmp_eta = fread(fid, nx * ny, 'double');
        if numel(tmp_eta) < nx * ny
            it = it_read - 1;
            break;
        end
        eta(it_read, :) = tmp_eta;
        fread(fid, 2, int_nbit);

        tmp_etax = fread(fid, nx * ny, 'double');
        if numel(tmp_etax) < nx * ny
            it = it_read - 1;
            break;
        end
        etax(it_read, :) = tmp_etax; %#ok<NASGU>
        fread(fid, 2, int_nbit);

        tmp_etay = fread(fid, nx * ny, 'double');
        if numel(tmp_etay) < nx * ny
            it = it_read - 1;
            break;
        end
        etay(it_read, :) = tmp_etay; %#ok<NASGU>
        fread(fid, 2, int_nbit);

        tmp_phi = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_phi) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        phi(it_read, :) = tmp_phi;
        fread(fid, 2, int_nbit);

        tmp_u = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_u) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        u(it_read, :) = tmp_u;
        fread(fid, 2, int_nbit);

        tmp_v = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_v) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        v(it_read, :) = tmp_v;
        fread(fid, 2, int_nbit);

        tmp_w = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_w) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        w(it_read, :) = tmp_w;
        fread(fid, 2, int_nbit);

        tmp_wz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_wz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        wz(it_read, :) = tmp_wz;
        fread(fid, 2, int_nbit);

        tmp_uz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_uz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        uz(it_read, :) = tmp_uz;
        fread(fid, 2, int_nbit);

        tmp_vz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_vz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        vz(it_read, :) = tmp_vz;
        fread(fid, 2, int_nbit);

        it = it_read;
    end

    if it <= 0
        error('No complete stored kinematics time step could be read from %s', file_path);
    end

    if it < nt
        eta = eta(1:it, :, :);
        phi = phi(1:it, :, :, :);
        u = u(1:it, :, :, :);
        v = v(1:it, :, :, :);
        w = w(1:it, :, :, :);
        uz = uz(1:it, :, :, :);
        vz = vz(1:it, :, :, :);
        wz = wz(1:it, :, :, :);
        t = t(1:it);
        nt = it;
    end

    if compute_derivatives
        alpha = 2;
        r = 2 * alpha + 1;
        c = build_stencil_even_local_package(alpha, 1);
        dt_matrix = spdiags(ones(nt, 1) * c(:, alpha + 1)', -alpha:alpha, nt, nt);
        for j = 1:alpha
            dt_matrix(j, :) = 0;
            dt_matrix(j, 1:r) = c(:, j)';
            dt_matrix(nt - j + 1, :) = 0;
            dt_matrix(nt - j + 1, nt - r + 1:nt) = c(:, r - j + 1)';
        end
        dt_matrix = dt_matrix / dt;

        etat_m = zeros(nt, size(eta, 2), size(eta, 3));
        etatt_m = zeros(nt, size(eta, 2), size(eta, 3));
        phit_m = zeros(size(phi));
        p_m = zeros(size(phi));
        ut_m = zeros(size(phi));

        for idy = 1:ny
            etat = zeros(nt, nx);
            etatt = zeros(nt, nx);
            phit = zeros(nt, nz, nx);
            p = zeros(nt, nz, nx);
            ut = zeros(nt, nz, nx);

            for ip = 1:nx
                eta_col = eta(:, ip, idy);
                etat(:, ip) = dt_matrix * eta_col;
                etatt(:, ip) = dt_matrix * etat(:, ip);

                for j = 1:nz
                    phi_col = phi(:, j, ip, idy);
                    w_col = w(:, j, ip, idy);
                    u_col = u(:, j, ip, idy);
                    uz_col = uz(:, j, ip, idy);

                    switch lower(phit_mode)
                        case 'uncorrected'
                            phit(:, j, ip) = dt_matrix * phi_col;
                        case 'sigma_corrected'
                            phit(:, j, ip) = dt_matrix * phi_col - w_col .* sigma(j) .* etat(:, ip);
                        otherwise
                            error('Unsupported phit_mode: %s', phit_mode);
                    end
                    p(:, j, ip) = -(phit(:, j, ip) + 0.5 * (u_col.^2 + v(:, j, ip, idy).^2 + w_col.^2));
                    ut(:, j, ip) = dt_matrix * u_col - uz_col .* sigma(j) .* etat(:, ip);
                end
            end

            etat_m(:, :, idy) = etat;
            etatt_m(:, :, idy) = etatt;
            phit_m(:, :, :, idy) = phit;
            p_m(:, :, :, idy) = p;
            ut_m(:, :, :, idy) = ut;
        end
    else
        etat_m = 0;
        etatt_m = 0;
        phit_m = 0;
        p_m = 0;
        ut_m = 0;
    end
end

function fx = build_stencil_even_local_package(alpha, der)
    rank = 2 * alpha + 1;
    fx = zeros(rank, rank);

    for ip = 1:alpha
        mat = zeros(rank, rank);
        row = 1;
        for m = -ip + 1:rank - ip
            for n = 1:rank
                mat(row, n) = m^(n - 1) / factorial(n - 1);
            end
            row = row + 1;
        end
        minv = inv(mat);
        fx(:, ip) = minv(der + 1, :).';
    end

    mat = zeros(rank, rank);
    row = 1;
    for m = -alpha:alpha
        for n = 1:rank
            mat(row, n) = m^(n - 1) / factorial(n - 1);
        end
        row = row + 1;
    end
    minv = inv(mat);
    fx(:, alpha + 1) = minv(der + 1, :).';

    if mod(der, 2) == 0
        for ip = 1:alpha
            fx(:, rank - ip + 1) = flipud(fx(:, ip));
        end
    else
        for ip = 1:alpha
            fx(:, rank - ip + 1) = -flipud(fx(:, ip));
        end
    end
end

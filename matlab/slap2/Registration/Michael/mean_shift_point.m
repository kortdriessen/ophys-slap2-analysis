function new_point = mean_shift_point(swc_point, im, radius, n_iter)
    if nargin < 4
        n_iter = 10;
    end

    cx = swc_point(2);
    cy = swc_point(1);
    cz = swc_point(3);

    % fprintf('Initial point: [%f, %f, %f]\n', cz, cy, cx);

    % Interpolate function from previous conversion
    interpolant = griddedInterpolant(im);

    count = 0;
    displacements = [];
    disp = Inf;
    while disp >= 0.5 && count < n_iter
        count = count + 1;

        % Sphere coordinates
        sc = sphere_coords(cx, cy, cz, radius);

        % Interpolated intensity values
        vals = interpolant(sc(:,[2 1 3]));

        % Compute weights
        w = vals - mean(vals);

        % Get positive weights and their indices
        w_idx = w > 0;
        if sum(w_idx) == 0
            break;
        end
        w_pos = w(w_idx);

        % Corresponding sphere coordinates
        c_pos = sc(w_idx, :);

        % Compute shifted coordinate
        c = sum(c_pos .* w_pos, 1) / sum(w_pos);

        % Displacement
        disp = norm(c - [cz, cy, cx]);
        displacements(end+1) = disp;

        % Update coordinates
        cz = c(3);
        cy = c(2);
        cx = c(1);
    end

    new_point = [cy, cx, cz];

    % fprintf('Adjusted: [%f, %f, %f]\n', cz, cy, cx);
end

function coords = sphere_coords(x, y, z, r)
    rr = r^2;
    coords = [];

    for dz = linspace(z - r, z + r, 2 * r + 1)
        for dy = linspace(y - r, y + r, 2 * r + 1)
            for dx = linspace(x - r, x + r, 2 * r + 1)
                dd = (dx - x)^2 + (dy - y)^2 + (dz - z)^2;
                if dd <= rr
                    coords = [coords; dx, dy, dz];
                end
            end
        end
    end
end


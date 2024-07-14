function guess = repairParameters(parameters, delete_position, gaussians, omega, chirp, first_color, second_color, gaussian_variation)
    % Receive parameters from a finished job, and reconstruct it to a
    % N_gaussians x 6 matrix

    for pos = flip(delete_position)
        parameters = [parameters(1:pos-1) 0 parameters(pos:end)];
    end

    if gaussian_variation
        % Handle different numbers of Gaussians for each color
        total_gaussians = first_color + second_color;
    else
        total_gaussians = gaussians;
    end

    parameters = reshape(parameters, total_gaussians, []);

    if (~chirp && size(parameters, 2) < 6)
        parameters = [parameters zeros(size(parameters, 1), 1)];
    end

    color_override = [];
    if gaussian_variation
        % Assign colors based on the number of Gaussians for each color
        for i = 1:first_color
            color_override = [color_override; omega(1)];
        end
        for i = 1:second_color
            color_override = [color_override; omega(2)];
        end
    else
        % Default behavior
        for color = omega
            for pulse = 1:size(parameters, 1) / size(omega, 2)
                color_override = [color_override; color];
            end
        end
    end

    if size(parameters, 2) < 6
        guess = [color_override parameters(:,:)];
    else
        guess = parameters;
    end
end

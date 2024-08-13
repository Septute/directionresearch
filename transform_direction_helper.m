% Helper function for transforming directions
function new_direction = transform_direction_helper(direction)
    if isnan(direction)
        new_direction = NaN;
    else
        direction = mod(direction, 360);
        if direction == 72
            new_direction = 90;
        elseif direction == 342 || direction == -18
            new_direction = -90;
        elseif direction == 0
            new_direction = -72;
        elseif direction == 180
            new_direction = 108;
        elseif direction >= 0 && direction < 72
            new_direction = -72 + (162 / 72) * direction;
        elseif direction > 72 && direction <= 180
            new_direction = 90 + (18 / 108) * (direction - 72);
        elseif direction > 180 && direction <= 360
            new_direction = 108 + (180 / 180) * (direction - 180);
        else
            new_direction = NaN;
        end
    end
end
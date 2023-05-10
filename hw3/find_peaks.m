function [] = find_peaks(data)
    fprintf("local maxima values \n");
    maxima_indices = find(diff(sign(diff(data))) < 0) + 1;
    local_maxima = data(maxima_indices);
    disp(local_maxima);

    fprintf("local minima values \n");
    minima_indices = find(diff(sign(diff(data))) > 0) + 1;
    local_minima = data(minima_indices);
    disp(local_minima);
end


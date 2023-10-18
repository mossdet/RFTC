function [xx, yy, zz, cc] = getScatterPlot3D_Vectors(x, y, z, c)
xx = []; yy = []; zz = []; cc = [];

    for xi = 1:length(x)
        for yi = 1:length(y)
            for zi = 1:length(z)
                xx = [xx x(xi)];
                yy = [yy y(yi)];
                zz = [zz z(zi)];
                cc = [cc c(zi)];
            end
        end
    end

end
function [x, y, dt] = simple_polygon2(numSides)

    if numSides < 3
        x = [];
        y = [];
        dt = DelaunayTri();
        return
    end

    oldState = warning('off', 'MATLAB:TriRep:PtsNotInTriWarnId');

    fudge = ceil(numSides/10);
    x = rand(numSides+fudge, 1);
    y = rand(numSides+fudge, 1);

    for i = 1:numSides+fudge
        for j = 1:numSides+fudge
            if i == j
                % do nothing
            else
                v = [x(j), y(j)] - [x(i), y(i)];
                if norm(v) < 0.7
                    x(j) = x(j) + 0.7*v(1)/(norm(v) + 1e-6);
                    y(j) = y(j) + 0.7*v(2)/(norm(v) + 1e-6);
                end
            end

        end
    end

    dt = DelaunayTri(x, y);
    boundaryEdges = freeBoundary(dt);
    numEdges = size(boundaryEdges, 1);

    while numEdges ~= numSides
        if numEdges > numSides
            triIndex = vertexAttachments(dt, boundaryEdges(:,1));
            triIndex = triIndex(randperm(numel(triIndex)));
            keep = (cellfun('size', triIndex, 2) ~= 1);
        end
        if (numEdges < numSides) || all(keep)
            triIndex = edgeAttachments(dt, boundaryEdges);
            triIndex = triIndex(randperm(numel(triIndex)));
            triPoints = dt([triIndex{:}], :);
            keep = all(ismember(triPoints, boundaryEdges(:,1)), 2);
        end
        if all(keep)
            warning('Couldnot achieve desired number of sides!');
            break
        end
        triPoints = dt.Triangulation;
        triPoints(triIndex{find(~keep, 1)}, :) = [];
        dt = TriRep(triPoints, x, y);
        boundaryEdges = freeBoundary(dt);
        numEdges = size(boundaryEdges, 1);
    end

    boundaryEdges = [boundaryEdges(:,1); boundaryEdges(1,1)];
    x = dt.X(boundaryEdges, 1);
    y = dt.X(boundaryEdges, 2);

    warning(oldState);

end
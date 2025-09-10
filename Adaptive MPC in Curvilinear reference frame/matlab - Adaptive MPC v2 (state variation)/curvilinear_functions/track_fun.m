function [PointAndTangent, x0] = track_fun(track_n, W_track)


%% Tracks description: spec = [length of curve, radius of curvature]
if track_n == 1
    x0 = [0; 0; pi/2];
    spec = [1, 0
            3*(pi/2), -3
            2, 0
            3*(pi/2), 3];
elseif track_n == 2
    x0 = [0; 0; 0];
    spec = [5, 0
            6*(pi/2), 6
            6*(pi/2), -6
            5, 0];
elseif track_n == 3
    x0 = [0; 0; 0];
    spec = [2, 0
            4*(pi/2), 4
            6, 0
            4*(pi/2), -4
            2, 0
            4*(pi/2), -4
            6, 0
            4*(pi/2), +4
            2, 0
            4*(pi/2), +4
            10, 0];
else
    fprintf("no track available");
end



%% PointAndTangent = [x, y, psi, s, length, curvature)
PointAndTangent = zeros(size(spec, 1), 6);
Centers = NaN(size(spec, 1), 2);


for i = 1:size(spec+1, 1)
    i;
    
    l = spec(i, 1);     % Length of the segments
    r = spec(i, 2);     % curvature radius

    if(spec(i, 2) == 0)     % straight path
        
        % fprintf("straight path");
        c = 0;
        
        %% Initial point
        if i == 1
            psi = x0(3);
            x = x0(1) + l*cos(psi);
            y = x0(2) + l*sin(psi);
            s = PointAndTangent(i, 4);
            
            Centers(i,:) = [(x+x0(1))/2, (y+x0(2))/2];
        else
            psi = PointAndTangent(i-1, 3);
            x = PointAndTangent(i-1, 1) + l*cos(psi);
            y = PointAndTangent(i-1, 2) + l*sin(psi);
            s = PointAndTangent(i-1, 4) + PointAndTangent(i-1, 5);
            
            Centers(i,:) = [(x+PointAndTangent(i-1, 1))/2, (y+PointAndTangent(i-1, 2))/2];
        end
        

        
        
    else                    % curve path
        % fprintf("curve path");
        c = 1/r;
        %% Curve direction
        if r >= 0
            direction = 1;
        else
            direction = -1;
        end
        
        
        %% Center of curvature
        if i == 1
            ang = x0(3);
            CenterX = x0(1) + abs(r) * cos(ang + direction * pi/2);
            CenterY = x0(2) + abs(r) * sin(ang + direction * pi/2);
            s = PointAndTangent(i, 4);
        else
            ang = PointAndTangent(i-1, 3);
            CenterX = PointAndTangent(i-1, 1) + abs(r)*cos(ang + direction*pi/2);
            CenterY = PointAndTangent(i-1, 2) + abs(r)*sin(ang + direction*pi/2);
            s = PointAndTangent(i-1, 4) + PointAndTangent(i-1, 5);
        end
        
        spanAng = l/abs(r);
        psi = wrap(ang + spanAng * sign(r));
        
        
        angleNormal = wrap(direction*pi/2 + ang);
        
        %% WIP
%         angle = -(pi - abs(angleNormal))*(sign(angleNormal));
        if angleNormal == 0
            angle = -(pi - abs(angleNormal));
        else
            angle = -(pi - abs(angleNormal))*(sign(angleNormal));
        end
        
        x = CenterX + abs(r)*cos(angle + direction*spanAng);
        y = CenterY + abs(r)*sin(angle + direction*spanAng);
        
        Centers(i,:) = [CenterX, CenterY];
        
        
    end
    
    
    %% PointAndTangent
    PointAndTangent(i, :) = [x, y, psi, s, l, c];
    

    
end



%% Plot reference
figure(1); hold on; grid on; axis equal;
title('Trajectory');
for i = 1:size(PointAndTangent,1)
    if i == 1
        plot([x0(1), PointAndTangent(1,1)], [x0(1), PointAndTangent(1,2)], '-k');
        plot([x0(1), PointAndTangent(1,1)], [x0(1), PointAndTangent(1,2)], '-k'); 
    elseif PointAndTangent(i,6) == 0
        plot([PointAndTangent(i-1,1), PointAndTangent(i,1)], [PointAndTangent(i-1,2), PointAndTangent(i,2)], '-k')
    else
        arc_center_and_edges([PointAndTangent(i-1,1); PointAndTangent(i-1,2)], [PointAndTangent(i,1); PointAndTangent(i,2)], Centers(i,:)', PointAndTangent(i,6), 0, 'k-')
        arc_center_and_edges([PointAndTangent(i-1,1); PointAndTangent(i-1,2)], [PointAndTangent(i,1); PointAndTangent(i,2)], Centers(i,:)', PointAndTangent(i,6), W_track, 'k--')
        arc_center_and_edges([PointAndTangent(i-1,1); PointAndTangent(i-1,2)], [PointAndTangent(i,1); PointAndTangent(i,2)], Centers(i,:)', PointAndTangent(i,6), -W_track, 'k--')
    end
end




return
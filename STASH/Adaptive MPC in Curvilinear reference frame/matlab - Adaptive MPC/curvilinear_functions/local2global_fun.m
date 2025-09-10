function [x, y] = local2global_fun(state, state0, PointAndTangent)


s = state(1);
ey = state(2);
epsi = state(3);

s0 = state0(1);
ey0 = state0(2);
epsi0 = state0(3);

TrackLength = PointAndTangent(end, 4) + PointAndTangent(end, 5);

while (s > TrackLength)
    s = s - TrackLength;
end


% Compute the segment in which system is evolving
i = find(s >= PointAndTangent(:,4) & s < PointAndTangent(:,4)+PointAndTangent(:,5), 1);

if (PointAndTangent(i, 6) == 0)     % segment is a straight line
    % extract the final point of the segment
    xf = PointAndTangent(i, 1);
    yf = PointAndTangent(i, 2);
    psi = PointAndTangent(i, 3);
    
    % extract the initial point of the segment
    if i == 1
        xs = s0;
        ys = ey0;
        psi_old = epsi0;
    else
        xs = PointAndTangent(i-1, 1);
        ys = PointAndTangent(i-1, 2);
        psi_old = PointAndTangent(i-1, 3);
    end
    
    
    % compute the segment length
    deltaL = PointAndTangent(i, 5);
    reltaL = s - PointAndTangent(i, 4);
    
    % linear combination
    x = (1 - reltaL / deltaL)*xs + reltaL/deltaL*xf + ey*cos(psi + pi/2);
    y = (1 - reltaL / deltaL)*ys + reltaL/deltaL*yf + ey*sin(psi + pi/2);

else                                    % segment is a curve
    r = 1/PointAndTangent(i, 6);        % Extract curvature
    ang = PointAndTangent(i - 1, 3);    % Extract angle of the tangent at the initial point (i-1)
    
    % Compute the center of the arc
    if r >= 0
        direction = 1;
    else
        direction = -1;
    end
    
    CenterX = PointAndTangent(i - 1, 1) + abs(r)*cos(ang + direction*pi/2);  % x coordinate center of circle
    CenterY = PointAndTangent(i - 1, 2) + abs(r)*sin(ang + direction*pi/2);  % y coordinate center of circle

    spanAng = (s - PointAndTangent(i, 4))/(pi*abs(r))*pi;
    angleNormal = wrap((direction*pi/2 + ang));
    angle = -(pi - abs(angleNormal))*(sign(angleNormal));
    
    x = CenterX + (abs(r) - direction*ey)*cos(angle + direction*spanAng);  % x coordinate of the last point of the segment
    y = CenterY + (abs(r) - direction*ey)*sin(angle + direction*spanAng);  % y coordinate of the last point of the segment

    
end


        
        

end
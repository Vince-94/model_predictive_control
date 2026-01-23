%% global2local
function [s, ey, epsi, CompletedFlag] = global2local_fun(state, state0, PointAndTangent, halfWidth)

s = 0;
ey = 0;
epsi = 0;
CompletedFlag = 0;

x = state(1);
y = state(2);
psi = state(3);

x0 = state0(1);
y0 = state0(2);
psi0 = state0(3);



for i = 1:size(PointAndTangent,1)
    i;
    
    if CompletedFlag == 1
        break;
    end
    
    %% Terminal points of segment
    % Initial segment point
    if i == 1
        xs = x0;
        ys = y0;
        psi_old = psi0;
    else
        xs = PointAndTangent(i-1, 1);
        ys = PointAndTangent(i-1, 2);
        psi_old = PointAndTangent(i-1, 3);
    end
    
    % Final segment point
    xf = PointAndTangent(i, 1);
    yf = PointAndTangent(i, 2);
    
    
    
    if PointAndTangent(i, 6) == 0
        %% epsi
        psi_unwrap = unwrap([psi], [], 1);
        epsi = psi_unwrap - psi_old;
        epsi = wrap(epsi);
        
        
        %% s, ey
        if (norm([xs, ys]-[x, y]) <= 10e-10)
            s = PointAndTangent(i, 4);
            ey = 0;
            CompletedFlag = 1;
        elseif (norm([xf, yf]-[x, y]) <= 10e-10)
            s = PointAndTangent(i, 4) + PointAndTangent(i, 5);
            ey = 0;
            CompletedFlag = 1;
        else
            if (abs(computeAngle([x,y], [xs,ys], [xf, yf])) <= pi/2+0.001 && abs(computeAngle([x,y], [xf,yf], [xs,ys])) <= pi/2+0.001)
                v1 = [x,y] - [xs, ys];
                angle = computeAngle([xf,yf], [xs,ys], [x,y]);
                s_local = norm(v1)*cos(angle);
                s = s_local + PointAndTangent(i, 4);
                ey = norm(v1)*sin(angle);

                if abs(ey) <= halfWidth
                    CompletedFlag = 1;
                end
                
            end

        end
            
        
    else
        
        %% Direction
        r = 1/PointAndTangent(i, 6);
%         r = PointAndTangent(i, 5)/PointAndTangent(i, 6);
        if r >= 0
            direction = 1;
        else
            direction = -1;
        end
        
        %% Arc center
        CenterX = xs + abs(r)*cos(psi_old + direction*pi/2);
        CenterY = ys + abs(r)*sin(psi_old + direction*pi/2); 
        
        
        %% s, ey, epsi
        if (norm([xs, ys] - [x, y]) <= 10e-10)
            s = PointAndTangent(i, 4);
            ey = 0;
            psi_unwrap = unwrap(psi, [], 1);
            epsi = psi_unwrap - psi_old;
            epsi = wrap(epsi);
            CompletedFlag = 1;
        elseif (norm([xf, yf] - [x, y]) <= 10e-10)
            s = PointAndTangent(i, 4) + PointAndTangent(i, 5);
            ey = 0;
            psi_unwrap = unwrap(psi, [], 1);
            epsi = psi_unwrap - PointAndTangent(i, 3);
            epsi = wrap(epsi);
            CompletedFlag = 1;
        else
            arc1 = PointAndTangent(i, 5)*PointAndTangent(i, 6);
            arc2 = computeAngle([xs,ys], [CenterX,CenterY], [x,y]);
            
            if (sign(arc1) == sign(arc2) && abs(arc1) >= abs(arc2))
                s_local = abs(arc2)*abs(r);
                s = s_local + PointAndTangent(i, 4);
                v = [x, y] - [CenterX, CenterY];
                ey = -sign(direction)*(norm(v) - abs(r));
                psi_unwrap = unwrap(psi, [], 1);
                epsi = psi_unwrap - (psi_old + arc2);
                epsi = wrap(epsi);
                
                if abs(ey) <= halfWidth
                    CompletedFlag = 1;
                end
            end
        end
        
    end
    

    
end


if CompletedFlag == 0
    s = 10000;
    ey = 10000;
    epsi = 10000;
    
    fprintf("Error!! POINT OUT OF THE TRACK!!!! \n")
%     dbstop if error 
end

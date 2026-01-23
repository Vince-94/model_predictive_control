function [angle] = computeAngle(point1, origin, point2)

v1 = point1 - origin;
v2 = point2 - origin;

dot = v1(1)*v2(1) + v1(2)*v2(2);
det = v1(1)*v2(2) - v1(2)*v2(1);
angle = atan2(det, dot);

end
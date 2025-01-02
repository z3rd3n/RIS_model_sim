function a_n = array_resp(theta,phi,x,y)
global k dis
a_n = exp(1i*k*dis*(x*sind(theta) + y*sind(phi)*cosd(theta) ));
end


utilities::utilities () {
  
  rot11 = cos(a) + (x * x) * (1 - cos(a));
  rot21 = z * sin(a) + x * y * (1 - cos(a));
  rot31 = y * sin(a) + x * z * (1 - cos(a));
  rot12 = x * y * (1 - cos(a)) - z * sin(a);
  rot22 = cos(a) + (y * y) * (1 - cos(a));
  rot32 = x * sin(a) + y * z * (1 - cos(a));
  rot13 = y * sin(a) + x * z * (1 - cos(a));
  rot23 = x * sin(a) + y * z * (1 - cos(a));
  rot33 = cos(a) + (z * x) * (1 - cos(a));
  
  rot23 = (-1) * rot23;
  rot31 = (-1) * rot31; 
  
}
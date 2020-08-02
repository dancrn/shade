vec4
fragment
(vec2 p)
{
  //your fancy shader code here
  float s = sin(uTime);
  return vec4(s*s, p.x, p.y, 1);
}

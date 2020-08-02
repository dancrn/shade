//based on the shader featured at:
//https://www.shadertoy.com/view/XdlSDs

vec4
fragment
(vec2 pt)
{
  vec2 p = (2.0*pt*uSize-uSize)/uSize.y;
  float tau = 3.1415926535*2.0;
  float a = atan(p.x,p.y);
  float r = length(p)*0.75;
  vec2 uv = vec2(a/tau,r);
  
  //get the color
  float xCol = (uv.x - (uTime / 3.0)) * 3.0;
  xCol = mod(xCol, 3.0);
  vec3 horColour = vec3(0.25, 0.25, 0.25);
  
  if (xCol < 1.0) {
    horColour.r += 1.0 - xCol;
    horColour.g += xCol;
  }
  else if (xCol < 2.0) {
    
    xCol -= 1.0;
    horColour.g += 1.0 - xCol;
    horColour.b += xCol;
  }
  else {
    xCol -= 2.0;
    horColour.b += 1.0 - xCol;
    horColour.r += xCol;
  }

  // draw color beam
  uv = (2.0 * uv) - 1.0;
  float beamWidth = (0.7+0.5*cos(uv.x*10.0*tau*0.15*clamp(floor(5.0 + 10.0*cos(uTime)), 0.0, 10.0))) * abs(1.0 / (30.0 * uv.y));
  return vec4((beamWidth * horColour), 1.0);
}

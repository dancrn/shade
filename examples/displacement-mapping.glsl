//based on the shader featured at:
//https://www.shadertoy.com/view/MtBSzd

//Displacement Mapping with ray tracing
#define PIXEL_SAMPLES     2
#define MAX_DISPLACEMENT  5.0
#define DISP_MAP_W        32
#define DISP_MAP_H        32
#define SHADOW
#define LIGHT_SAMPLES     2
#define GAMMA             2.2
#define SIN_WAVE

float displacementMultiplier;

//used macros and constants
#define HALF_PI             1.5707963
#define PI                  3.1415926
#define TWO_PI              6.2831852
#define FOUR_PI             12.566370
#define INV_PI              0.3183099
#define INV_TWO_PI          0.1591549
#define INV_FOUR_PI         0.0795775
#define EPSILON             0.0001
#define IN_RANGE(x, a, b)   (((x) > (a)) && ((x) < (b)))
#define EQUAL_FLT(a, b, q)  (((a)>((b)-(q))) && ((a)<((b)+(q))))
#define IS_ZERO(a)          EQUAL_FLT(a, 0.0, EPSILON)

//********************************************
// random number generator **********
// taken from iq :)
float seed;  //seed initialized in main
float rnd() { return fract(sin(seed++)*43758.5453123); }
//***********************************

// Data structures ****************** 
struct AABB { vec3 min_; vec3 max_; };
struct Triangle { vec3 v0; vec3 v1; vec3 v2; };
struct Range { float min_; float max_; };
struct Ray { vec3 origin; vec3 dir; };
struct Camera { mat3 rotate; vec3 pos; float fovV; };
//***********************************
AABB displacementVolume;// = AABB(vec3(-5.0,-5.0, 0.0), vec3(5.0,5.0, MAX_DISPLACEMENT));

// ************************  INTERSECTION FUNCTIONS **************************
bool rayAABBIntersection( in Ray ray, AABB aabb, out float t_enter, out float t_exit ) {
    vec3 OMIN = ( aabb.min_ - ray.origin ) / ray.dir;
    vec3 OMAX = ( aabb.max_ - ray.origin ) / ray.dir;
    vec3 MAX = max ( OMAX, OMIN );
    vec3 MIN = min ( OMAX, OMIN );
    t_exit = min ( MAX.x, min ( MAX.y, MAX.z ) );
    t_enter = max ( MIN.x, max ( MIN.y, MIN.z ) );
    
    if ( t_exit <= t_enter )
        return false;
    
    return ( t_exit > t_enter );
}

vec3 getTriangleNormal(vec3 v0, vec3 v1, vec3 v2) {
    return normalize(cross(v1-v0,v2-v0));
}

bool rayIntersectsTriangle(vec3 p, vec3 d, vec3 v0, vec3 v1, vec3 v2, out float t){
  vec3 e1,e2,h,s,q;
  float a,f,u,v;
  e1 = v1-v0;
  e2 = v2-v0;

  h = cross(d,e2);
  a = dot(e1,h);

  if (a > -0.00001 && a < 0.00001)
    return false;

  f = 1.0 / a;
  s = p-v0;
  u = f * dot(s,h);

  if (u < 0.0 || u > 1.0)
    return false;

  q = cross(s,e1);
  v = f * dot(d,q);

  if (v < 0.0 || u + v > 1.0)
    return false;

  // at this stage we can compute t to find out where
  // the intersection point is on the line
  t = f * dot(e2,q);

  //uv = vec2(u, v);
  
  if (t > 0.00001) // ray intersection
    return true;

  // this means that there is a line intersection
  // but not a ray intersection
  return false;
}

bool rayQuadIntersect( in Ray ray, vec3 v0, vec3 v1, vec3 v2, vec3 v3, out float t, out vec3 n ) {
    float tcurrent;
    t = 10e+10;
  
    //first triangle
    if ( rayIntersectsTriangle( ray.origin, ray.dir, v0, v1, v2, tcurrent) ) {
        t = tcurrent;
        n = getTriangleNormal(v0, v1, v2);
    }
    
    //second triangle
    if ( rayIntersectsTriangle( ray.origin, ray.dir, v0, v2, v3, tcurrent) ) {
        if(tcurrent<t) {
            t = tcurrent;
          n = getTriangleNormal(v0, v2, v3);
        }
    }
    
    return (t < 10e+10);
}
// ***************************************************************************

vec3 uniformDirectionWithinCone( in vec3 d, in float phi, in float sina, in float cosa ) {    
  vec3 w = normalize(d);
    vec3 u = normalize(cross(w.yzx, w));
    vec3 v = cross(w, u);
  return (u*cos(phi) + v*sin(phi)) * sina + w * cosa;
}

///////////////////////////////////////////////////////////////////////
void initCamera( in vec3 pos, in vec3 frontDir, in vec3 upDir, in float fovV, out Camera dst ) {
  vec3 back = normalize( -frontDir );
  vec3 right = normalize( cross( upDir, back ) );
  vec3 up = cross( back, right );
    dst.rotate[0] = right;
    dst.rotate[1] = up;
    dst.rotate[2] = back;
    dst.fovV = fovV;
    dst.pos = pos;
}

Ray genRay( in Camera camera, in vec2 pixel ) {
  vec2 uResolution = vec2(1024, 768);
  vec2 iPlaneSize=2.*tan(0.5*camera.fovV)*vec2(uResolution.x/uResolution.y,1.);
  vec2 ixy=(pixel/uResolution.xy - 0.5)*iPlaneSize;
    
    Ray ray;
    ray.origin = camera.pos;
  ray.dir = camera.rotate*normalize(vec3(ixy.x,ixy.y,-1.0));

  return ray;
}

#define PIXEL2UV(x,y) vec2( float(x)/float(DISP_MAP_W), float(y)/float(DISP_MAP_H) )
float getDisplacement(vec2 uv){
    
    return (sin(uTime-length(vec2(0.5,0.5) - uv)*20.0)+1.0)/2.0;
}

void getPixelDisplacements( vec2 uv, out float d1, out float d2, out float d3, out float d4 ) {
    int px = int(uv.x*float(DISP_MAP_W));
    int py = int(uv.y*float(DISP_MAP_H));

    d1 = getDisplacement( PIXEL2UV(px  ,py  ) );
    d2 = getDisplacement( PIXEL2UV(px+1,py  ) );
    d3 = getDisplacement( PIXEL2UV(px+1,py+1) );
    d4 = getDisplacement( PIXEL2UV(px  ,py+1) );
}

#define UVH2POS(aabb,uvd,pos) { vec3 aabbDim = aabb.max_ - aabb.min_; pos = aabb.min_ + uvd*aabbDim; }
#define POS2UVH(aabb,pos,uvd) { vec3 aabbDim = aabb.max_ - aabb.min_; uvd = (pos-aabb.min_)/aabbDim; }
bool processVoxel(Ray ray, vec2 p, out float t, out vec3 normal) {
    //Lookup displacement values
    vec2 uv = p/vec2(DISP_MAP_W,DISP_MAP_H);
            
    float disp[4];
    getPixelDisplacements( uv, disp[0], disp[1], disp[2], disp[3] );
    
     //calculate displaced vertices
     int px = int(p.x);
    int py = int(p.y);

    vec3 corner_uv[4];
    corner_uv[0] = vec3(PIXEL2UV(px  ,py  ), disp[0]);
    corner_uv[1] = vec3(PIXEL2UV(px+1,py  ), disp[1]);
    corner_uv[2] = vec3(PIXEL2UV(px+1,py+1), disp[2]);
    corner_uv[3] = vec3(PIXEL2UV(px  ,py+1), disp[3]);

    vec3 vertices[4];
    UVH2POS(displacementVolume,corner_uv[0],vertices[0]);
    UVH2POS(displacementVolume,corner_uv[1],vertices[1]);
    UVH2POS(displacementVolume,corner_uv[2],vertices[2]);
    UVH2POS(displacementVolume,corner_uv[3],vertices[3]);

    float hitDist;
    vec3 hitN;
    if( rayQuadIntersect( ray, vertices[0], vertices[1], vertices[2], vertices[3], hitDist, hitN ) ) {
        t = hitDist;
        normal = hitN;
        return true;
    }
    
    return false;
}

bool processVoxelsOnLine(Ray ray, vec3 p0, vec3 p1, out float t, out vec3 normal) {
  int gx0idx = int(floor(p0.x));
  int gy0idx = int(floor(p0.y));
  int gz0idx = int(floor(p0.z));
  
  int gx1idx = int(floor(p1.x));
  int gy1idx = int(floor(p1.y));
  int gz1idx = int(floor(p1.z));
    
    int sx = (gx1idx > gx0idx) ? 1 : (gx1idx < gx0idx) ? -1 : 0;
    int sy = (gy1idx > gy0idx) ? 1 : (gy1idx < gy0idx) ? -1 : 0;
    int sz = (gz1idx > gz0idx) ? 1 : (gz1idx < gz0idx) ? -1 : 0;
        
  int gx = gx0idx;
  int gy = gy0idx;
  int gz = gz0idx;
  
  //Planes for each axis that we will next cross
    int gxp = gx0idx + (gx1idx > gx0idx ? 1 : 0);
    int gyp = gy0idx + (gy1idx > gy0idx ? 1 : 0);
    int gzp = gz0idx + (gz1idx > gz0idx ? 1 : 0);
    
  //Only used for multiplying up the error margins
  float vx = (p1.x == p0.x)? 1.0 : p1.x - p0.x;
  float vy = (p1.y == p0.y)? 1.0 : p1.y - p0.y;
  float vz = (p1.z == p0.z)? 1.0 : p1.z - p0.z;
  
    //Error is normalized to vx * vy * vz so we only have to multiply up
    float vxvy = vx * vy;
    float vxvz = vx * vz;
    float vyvz = vy * vz;
  
  //Error from the next plane accumulators, scaled up by vx*vy*vz
  float errx = (float(gxp) - p0.x) * vyvz;
  float erry = (float(gyp) - p0.y) * vxvz;
  float errz = (float(gzp) - p0.z) * vxvy;
  
  float derrx = float(sx) * vyvz;
  float derry = float(sy) * vxvz;
  float derrz = float(sz) * vxvy;
    
    for(int i=0; i<(DISP_MAP_W+DISP_MAP_H); i++) {
        if( processVoxel(ray, vec2(float(gx),float(gy)), t, normal) ) {
            return true;
        }
    
    if ((gx == gx1idx) && (gy == gy1idx) && gz == gz1idx)
            break;

        //Which plane do we cross first?
    float xr = abs(errx);
    float yr = abs(erry);
    float zr = abs(errz);
        
        //console.log("err",errx,erry,errz);
    if ((sx != 0) && (sy == 0 || xr < yr) && (sz == 0 || xr < zr)) {
      gx += sx;
      errx += derrx;
    } else if (sy != 0 && (sz == 0 || yr < zr)) {
      gy += sy;
      erry += derry;
    } else if (sz != 0) {
      gz += sz;
      errz += derrz;
    }
  }
    return false;
}

bool rayIntersectsDisplacement( in Ray ray, out float t, out vec3 normal, out int iter ) {
    iter = 0;
    float t1, t2;
    if( rayAABBIntersection( ray, displacementVolume, t1, t2) ) {
         vec3 hitpos1 = ray.origin + ray.dir*(t1+EPSILON);  //volume entry point
        vec3 hitpos2 = ray.origin + ray.dir*(t2-EPSILON);  //volume exit point
        
        //Convert position to parametric coordinates
        vec3 uvd1, uvd2;
        POS2UVH(displacementVolume,hitpos1,uvd1);
        POS2UVH(displacementVolume,hitpos2,uvd2);
        
        //pixel coordinates of projected entry and exit point
        vec3 p0 = uvd1*vec3(float(DISP_MAP_W),float(DISP_MAP_H),1.0);
        vec3 p1 = uvd2*vec3(float(DISP_MAP_W),float(DISP_MAP_H),1.0);
        
        if( processVoxelsOnLine(ray, p0, p1, t, normal) ) {
            iter = 1;
            return true;
        }
    }
    
    return false;
}

vec3 getColor(Ray ray) {
    float t;
    vec3 n;
    int iter;
    if( rayIntersectsDisplacement( ray, t, n, iter ) ) {
        vec3 p = ray.origin + ray.dir*t;
        vec3 normal = dot(ray.dir,n)<0.0?n:-n;
        
        vec3 lc = normalize(vec3(0.4,0.5,0.3));
        float a_max = PI/30.0;
        
        
    vec3 Li = vec3(0.0);
        for(int i=0; i<LIGHT_SAMPLES; i++) {
            float a = rnd()*a_max;
            float sina = sin(a);
          float cosa = sqrt(1.0-sina*sina);
            vec3 l = uniformDirectionWithinCone( lc, rnd()*TWO_PI, sina, cosa );

            float dotNL = dot(normal,l);
            if(dotNL < 0.0) {
                return abs(dotNL)*vec3(0.01);//fake indirect
            }
            
            Ray shadowRay = Ray( p + normal*EPSILON, l );

#ifdef SHADOW
            if(rayIntersectsDisplacement( shadowRay, t, n, iter )) {
              Li += vec3(0.003);//fake indirect
            }else
#endif
            {
              Li += max(0.0, dotNL)*vec3(1.0);
            }
            
      }
        
        return Li*(1.0/float(LIGHT_SAMPLES));
    }
    return vec3(0.0);
}

vec4
fragment
(vec2 p)
{
  vec2 uResolution = vec2(1024, 768);
  vec4 uMouse = vec4(0, 0, 0, 0);
  seed = uResolution.y * p.x / uResolution.x + p.y / uResolution.y;
    
  float sinTime = sin(uTime*0.2);
  vec3 cameraPos = vec3( 6.0 + sin(uTime*0.4), -9.0 + sin(uTime*0.3), 7.0 + sin(uTime*0.4)*5.0 );
  vec3 cameraTarget = vec3( 0.0, 0.0, 0.0 );
    
  displacementVolume.min_ = vec3(-5.0,-5.0, 0.0);
    displacementMultiplier = (uMouse.w!=0.0)?min(1.0,uMouse.y/uResolution.y):0.4;
    displacementVolume.max_ = vec3(5.0,5.0, displacementMultiplier*MAX_DISPLACEMENT);
    
    Camera camera;
   
    initCamera( cameraPos, cameraTarget - cameraPos, vec3( 0.0, 0.0, 1.0 ), radians(40.0), camera );
  
    Ray ray;
  vec3 accumulatedColor = vec3( 0.0 );
  for(int si=0; si<PIXEL_SAMPLES; ++si ){
        vec2 screenCoord = p*uResolution + vec2( (1.0/float(PIXEL_SAMPLES))*(float(si)+rnd()), rnd() );
        ray = genRay( camera, screenCoord );
        
        accumulatedColor += getColor( ray );
  }
  
  //devide to sample count
  accumulatedColor = accumulatedColor*(1.0/float(PIXEL_SAMPLES));
  
  //gamma correction
    accumulatedColor = pow( accumulatedColor, vec3( 1.0 / GAMMA ) );
    
  return vec4( accumulatedColor,1.0 );
}

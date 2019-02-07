#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
in vec4 fs_LightVec;
out vec4 out_Col;

vec3 u_LightColor = vec3(1.0);
float EPSILON = 0.1;

vec3 rayCast() {
  float sx = fs_Pos.x;
  float sy = fs_Pos.y;
  float len = length(u_Ref - u_Eye);
  vec3 look = normalize(u_Ref - u_Eye);
  vec3 right = normalize(cross(look, u_Up));
  vec3 up = cross(right, look);
  float tan_fovy = tan(45.0 / 2.0);
  float aspect = u_Dimensions.x / u_Dimensions.y;
  vec3 V = up * len * tan_fovy;
  vec3 H = right * len * aspect * tan_fovy;

  vec3 p = u_Ref + sx * H + sy * V;
  vec3 direct = normalize(p - u_Eye);
  return direct;
}

float sdPlane( vec3 p ) {
  return p.y;
}

float sdSphere( vec3 p, float s ) {
  return length(p)-s;
}

float sdBox( vec3 p, vec3 b ) {
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

vec2 opU( vec2 d1, vec2 d2 ) {
  return (d1.x<d2.x) ? d1 : d2;
}


vec2 map(vec3 pos) {
  vec2 res = opU( vec2( sdPlane(pos), 1.0 ),
  vec2( sdSphere(pos-vec3( 0.0,0.25, 0.0), 0.25 ), 46.9 ) );
  res = opU(res, vec2( sdBox(pos-vec3( 1.0,0.25, 0.0), vec3(0.25) ), 3.0 ) );
  return res;
}

vec2 castRay( in vec3 ro, in vec3 rd )
{
  float tmin = 1.0;
  float tmax = 20.0;

  #if 1
  // bounding volume
  float tp1 = (0.0-ro.y)/rd.y; if( tp1>0.0 ) tmax = min( tmax, tp1 );
  float tp2 = (1.6-ro.y)/rd.y; if( tp2>0.0 ) { if( ro.y>1.6 ) tmin = max( tmin, tp2 );
    else           tmax = min( tmax, tp2 ); }
    #endif

    float t = tmin;
    float m = -1.0;
    for( int i=0; i<64; i++ )
    {
      float precis = 0.0004*t;
      vec2 res = map( ro+rd*t );
      if( res.x<precis || t>tmax ) break;
      t += res.x;
      m = res.y;
    }

    if( t>tmax ) m=-1.0;
    return vec2( t, m );
  }

  vec3 sphereNormal(vec3 p, float rad) {
    return p / rad;
  }

  vec3 estimateNormal(vec3 pos) {
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*map( pos + e.xyy ).x +
    e.yyx*map( pos + e.yyx ).x +
    e.yxy*map( pos + e.yxy ).x +
    e.xxx*map( pos + e.xxx ).x );
  }

  vec3 lambert(vec3 normal, vec3 direction, vec3 color) {
    float diffuseTerm = dot(normalize(normal), normalize(direction));
    diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);
    float ambientTerm = 0.2;

    float lightIntensity = diffuseTerm + ambientTerm;
    return clamp(vec3(color.rgb * lightIntensity * u_LightColor), 0.0, 1.0);
  }


  void main() {
    float t1 = 0.01;
    vec4 color = vec4(1.0);
    vec3 p = u_Eye + t1 * rayCast();
    color = vec4(0.5 * (rayCast() + vec3(1.0, 1.0, 1.0)), 1.0);
    vec3 ro = u_Eye;
    vec3 rd = rayCast();
    vec3 col = vec3(0.7, 0.9, 1.0) + rd.y * 0.8;
    vec2 res = castRay(ro,rd);
    float t = res.x;
    float m = res.y;
    if( m>-0.5 )
    {
      vec3 pos = ro + t*rd;
      vec3 nor = estimateNormal(pos);
      col = lambert(nor, vec3(fs_LightVec), vec3(0.2, 0.3, 0.2));
    }

    // if (sdFloor(p) < 0.1) {
    //   vec3 col = lambert(estimateNormal(p), vec3(fs_LightVec), vec3(0.2, 1.0, 0.5));
    //   color = vec4(col, 1.0);
    //   break;
    // }

    out_Col = vec4(col, 1.0);
  }

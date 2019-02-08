#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform float u_Size;
uniform float u_Color;

in vec2 fs_Pos;
in vec4 fs_LightVec;
out vec4 out_Col;

vec3 u_LightColor = vec3(1.0);
float EPSILON = 0.1;

 struct BoundingBox
 {
    vec3 max;
    vec3 min;
 };

// noise functions
float random1( vec2 p , vec2 seed) {
  return fract(sin(dot(p + seed, vec2(127.1, 311.7))) * 43758.5453);
}

// BoundingBox Intersection test
// adapted from
//https://tavianator.com/fast-branchless-raybounding-box-intersections/
// and CIS
float boxIntersection(BoundingBox b, vec3 rayDir, vec3 rayOri) {
    float tx1 = (b.min.x - rayOri.x)* rayDir.x;
    float tx2 = (b.max.x - rayOri.x)* rayDir.x;

    float tmin = min(tx1, tx2);
    float tmax = max(tx1, tx2);

    float ty1 = (b.min.y - rayOri.y)* rayDir.y;
    float ty2 = (b.max.y - rayOri.y)* rayDir.y;

    tmin = max(tmin, min(ty1, ty2));
    tmax = min(tmax, max(ty1, ty2));

    float tz1 = (b.min.z - rayOri.z)* rayDir.z;
    float tz2 = (b.max.z - rayOri.z)* rayDir.z;

    tmin = max(tmin, min(tz1, tz2));
    tmax = min(tmax, max(tz1, tz2));

    if (tmax > tmin) {
      return tmax;
    }

    else {
      return 1000.0;
    }
}

float boundingVolumeHierachy(vec3 point, vec3 dir, BoundingBox boxes[2]) {
  float t = 1000.;
  for (int i = 0; i < boxes.length(); i++) {
		float intersect = boxIntersection(boxes[i], point, dir);
		if (t > intersect) {
			t = intersect;
		}
	}
  return t;
}

float interpNoise2D(float x, float y) {
  float intX = floor(x);
  float fractX = fract(x);
  float intY = floor(y);
  float fractY = fract(y);

  float v1 = random1(vec2(intX, intY), vec2(311.7, 127.1));
  float v2 = random1(vec2(intX + 1.0f, intY), vec2(311.7, 127.1));
  float v3 = random1(vec2(intX, intY + 1.0f), vec2(311.7, 127.1));
  float v4 = random1(vec2(intX + 1.0, intY + 1.0), vec2(311.7, 127.1));

  float i1 = mix(v1, v2, fractX);
  float i2 = mix(v3, v4, fractX);

  return mix(i1, i2, fractY);
}

float generateColor(float x, float y) {
  float total = 0.0;
  float persistence = 0.5f;
  float octaves = 10.0;

  for (float i = 0.0; i < octaves; i = i + 1.0) {
    float freq = pow(2.0f, i);
    float amp = pow(persistence, i);
    total += (1.0 / freq) * interpNoise2D(x * freq, y * freq);
  }
  return total;
}

// Raycasting
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

// primitives courtesy of http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

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

float sdCappedCone( in vec3 p, in float h, in float r1, in float r2 ) {
  vec2 q = vec2( length(p.xz), p.y );
  vec2 k1 = vec2(r2,h);
  vec2 k2 = vec2(r2-r1,2.0*h);
  vec2 ca = vec2(q.x-min(q.x,(q.y < 0.0)?r1:r2), abs(q.y)-h);
  vec2 cb = q - k1 + k2*clamp(dot(k1-q,k2)/ dot(k2, k2), 0.0, 1.0 );
  float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
  return s*sqrt( min(dot(ca, ca), dot(cb, cb)));
}

float sdTorus( vec3 p, vec2 t) {
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
  vec2 d = vec2( length(p.xz)-2.0*ra+rb, abs(p.y) - h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}

float sdHexPrism(vec3 p, vec2 h)
{
  const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
  p = abs(p);
  p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
  vec2 d = vec2(
    length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
    p.z-h.y );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
  }

  // operations courtesy of http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
  float opUnion( float d1, float d2 ) {
    return min(d1,d2);
  }

  vec2 opU( vec2 d1, vec2 d2 ) {
    return (d1.x<d2.x) ? d1 : d2;
  }

  float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h);
  }

  float opIntersection( float d1, float d2 ) { return max(d1,d2);}

  float opSubtraction( float d1, float d2 ) { return max(-d1,d2); }

  vec3 opTwist(vec3 p)
  {
    float  c = cos(10.0*p.y+10.0);
    float  s = sin(10.0*p.y+10.0);
    mat2   m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
  }

  float opOnion( in float sdf, in float thickness )
  {
    return abs(sdf)-thickness;
  }

  // ease in ease out
  float ease_in (float t) {
    return t * t;
  }

  float ease_in_ease_out (float t) {
    if (t < 0.5) {
      return ease_in(t * 2.0) / 2.0;
    }
    else {
      return 1.0 - ease_in((1.0 - t) * 2.0) / 2.0;
    }
  }

  // adding all the primitives into one map

  vec2 trashCan(vec3 pos) {
    vec2 res = opU(vec2(sdPlane(pos), 0.5),
    vec2(sdCappedCone(pos, 2.2, 0.1, 0.8), 25.0));
    res.x = opSmoothUnion(res.x, sdTorus(pos-vec3(0.0, 2.2, 0.0), vec2(0.8, 0.1)), .05);
    return res;
  }

  vec2 trash1(vec3 pos) {
    pos.y = pos.y - abs(-sin(u_Time/10.0));
    float t1 = 0.5*sdTorus(opTwist(pos-vec3(1.0, 0.25, 0.2)),vec2(0.20 + 0.05 * u_Size,0.05));
    float subtract = opSubtraction(t1, sdSphere(pos-vec3(1.0, 0.25, 0.2),
    0.20 + 0.05 * u_Size));
    return vec2(subtract, 0.5);
  }

  vec2 trash2(vec3 pos) {
    pos.y = pos.y - ease_in_ease_out(-sin(u_Time/15.0));
    float t2 = opIntersection((sdSphere(pos-vec3(0.7, 0.15, 0.1), 0.15 + 0.03 * u_Size)),
    sdBox(pos-vec3(0.7, 0.15, 0.1), vec3(0.10 + 0.03 * u_Size)));
    return vec2(t2, 0.5);
  }

  vec2 castRay_TrashCan( in vec3 ro, in vec3 rd, float t1)
  {
    float tmin = t1;
    float tmax = 20.0;

    float t = tmin;
    float m = -1.0;
    for( int i=0; i<64; i++ )
    {
      float precis = 0.0004*t;
      vec2 res = trashCan(ro + rd * t );
      if( res.x<precis || t>tmax ) break;
      t += res.x;
      m = res.y;
    }

    if( t>tmax ) m=-1.0;
    return vec2( t, m );
  }

  vec2 castRay_Trash1( in vec3 ro, in vec3 rd, float t1)
  {
    float tmin = t1;
    float tmax = 20.0;

    float t = tmin;
    float m = -1.0;
    for( int i=0; i<64; i++ )
    {
      float precis = 0.0004*t;
      vec2 res = trash1(ro + rd * t );
      if( res.x<precis || t>tmax ) break;
      t += res.x;
      m = res.y;
    }

    if( t>tmax ) m=-1.0;
    return vec2( t, m );
  }

  vec2 castRay_Trash2( in vec3 ro, in vec3 rd, float t1)
  {
    float tmin = t1;
    float tmax = 20.0;

    float t = tmin;
    float m = -1.0;
    for( int i=0; i<64; i++ )
    {
      float precis = 0.0004*t;
      vec2 res = trash2(ro + rd * t );
      if( res.x<precis || t>tmax ) break;
      t += res.x;
      m = res.y;
    }

    if( t>tmax ) m=-1.0;
    return vec2( t, m );
  }


  // calculating normals

  vec3 estimateNormal_TrashCan(vec3 pos) {
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*trashCan( pos + e.xyy ).x +
    e.yyx*trashCan( pos + e.yyx ).x +
    e.yxy*trashCan( pos + e.yxy ).x +
    e.xxx*trashCan( pos + e.xxx ).x);
  }

  vec3 estimateNormal_Trash1(vec3 pos) {
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*trash1( pos + e.xyy ).x +
    e.yyx*trash1( pos + e.yyx ).x +
    e.yxy*trash1( pos + e.yxy ).x +
    e.xxx*trash1( pos + e.xxx ).x);
  }

  vec3 estimateNormal_Trash2(vec3 pos) {
    vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
    return normalize( e.xyy*trash2( pos + e.xyy ).x +
    e.yyx*trash2( pos + e.yyx ).x +
    e.yxy*trash2( pos + e.yxy ).x +
    e.xxx*trash2( pos + e.xxx ).x);
  }

  // calculate lambert shading

  vec3 lambert(vec3 normal, vec3 direction, vec3 color) {
    float diffuseTerm = dot(normalize(normal), normalize(direction));
    diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);
    float ambientTerm = 0.2;
    float lightIntensity = diffuseTerm + ambientTerm;
    return clamp(vec3(color.rgb * lightIntensity * u_LightColor), 0.0, 1.0);
  }


  void main() {
    // trash
    BoundingBox boxes[2];

    BoundingBox trash;
    trash.min = vec3(-2.0, -1.0, -10.0);
    trash.max = vec3(5.0, 10.0, 10.0);

    BoundingBox ground;
    ground.min = vec3(-30.0, -10.0, -30.0);
    ground.max = vec3(30.0, 10.0, 30.0);

    boxes[0] = trash;
    boxes[1] = ground;

    vec4 color = vec4(1.0);
    float t_back = 0.01;
    vec3 p = u_Eye + t_back * rayCast();
    color = vec4(0.5 * (rayCast() + vec3(1.0, 1.0, 1.0)), 1.0)
    + vec4(vec3(0.10, 0.0, 0.0) * u_Color, 1.0);
    vec3 ro = u_Eye;
    vec3 rd = rayCast();
    float t1 = boundingVolumeHierachy(ro, rd, boxes);

    vec2 res1 = castRay_TrashCan(ro,rd, t1);
    float t = res1.x;
    float m = res1.y;

    if(m > -0.2)
    {
      vec3 pos = ro + t*rd;
      vec3 nor = estimateNormal_TrashCan(pos);
      color = vec4(lambert(nor, vec3(u_Eye), vec3(0.76)), 1.0);
      if (m < 1.0) {
        vec3 noiseColor = vec3(255. / 255., 247. / 255., 234. / 255.)
        * generateColor(pos.x, pos.z);
        color = vec4(lambert(nor, vec3(u_Eye), noiseColor), 1.0);
      }
    }

    vec2 res2 = castRay_Trash1(ro,rd, t1);
    float t2 = res2.x;
    float m2 = res2.y;
    if(m2 > -0.2)
    {
      vec3 pos = ro + t2*rd;
      vec3 nor = estimateNormal_Trash1(pos);
      color = vec4(lambert(nor, vec3(u_Eye), vec3(206. / 255., 145./255., 120./255.)), 1.0);
    }

    vec2 res3 = castRay_Trash2(ro,rd, t1);
    float t3 = res3.x;
    float m3 = res3.y;
    if(m3 > -0.2)
    {
      vec3 pos = ro + t3 * rd;
      vec3 nor = estimateNormal_Trash2(pos);
      color = vec4(lambert(nor, vec3(u_Eye), vec3(199./255., 219./255., 255./255.)), 1.0);
    }

    out_Col = color;
  }

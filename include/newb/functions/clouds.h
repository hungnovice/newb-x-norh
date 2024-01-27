#ifndef CLOUDS_H
#define CLOUDS_H

#include "noise.h"

// simple clouds noise
float cloudNoise2D(vec2 p, highp float t, float rain) {
  t *= NL_CLOUD1_SPEED;
  p += t;
  p.x += sin(p.y*0.4 + t);

  vec2 p0 = floor(p);
  vec2 u = p-p0;
  u *= u*(3.0-2.0*u);
  vec2 v = 1.0-u;

  // rain transition
  vec2 d = vec2(0.09+0.5*rain,0.089+0.5*rain*rain);

  return v.y*(randt(p0,d)*v.x + randt(p0+vec2(1.0,0.0),d)*u.x) +
         u.y*(randt(p0+vec2(0.0,1.0),d)*v.x + randt(p0+vec2(1.0,1.0),d)*u.x);
}

// simple clouds
vec4 render_clouds_simple(vec3 pos, highp float t, float rain, vec3 zenith_col, vec3 horizon_col, vec3 fog_col) {
  pos.xz *= NL_CLOUD1_SCALE;

  float cloudAlpha = cloudNoise2D(pos.xz, t, rain);
  float cloudShadow = cloudNoise2D(pos.xz*0.91, t, rain);

  vec4 color = vec4(0.02,0.04,0.05,cloudAlpha);

  color.rgb += fog_col;
  color.rgb *= 1.0-0.5*cloudShadow*step(0.0,pos.y);

  color.rgb += zenith_col*0.7;
  color.rgb *= 1.0 - 0.4*rain;

  return color;
}

// rounded clouds

//   apply bevel with radius r at at corner (1.0)
float bevel(float x, float r) {
	 //return smoothstep(0.5,1.5,x);
  float y = max(x-r,0.0);
  return r+sqrt(1.0-2.0*r+r*r-y*y);
}

// 2d noise
float noise(vec2 p){
  vec2 p0 = floor(p);
  vec2 u = p-p0;

  u *= u*(3.0-2.0*u);
  vec2 v = 1.0 - u;

  float c1 = rand(p0);
  float c2 = rand(p0+vec2(1.0,0.0));
  float c3 = rand(p0+vec2(0.0,1.0));
  float c4 = rand(p0+vec2(1.0,1.0));

  float n = v.y*mix(c1,c2,u.x) + u.y*(c3*v.x+c4*u.x);
  return n;
}

float cloud_df(vec3 pos, float rain) {
  pos.x += 0.17*noise(6.0*pos.xz);
  pos.y += 0.1*noise(0.0*pos.xz);
  pos.z += 0.17*noise(6.0*pos.xz);
  vec2 p0 = floor(pos.xz);
  vec2 u = pos.xz - p0;
  //u = smoothstep(0.99*NL_CLOUD2_SHAPE,1.0,u);
  u = max((u-NL_CLOUD2_SHAPE)/(1.0-NL_CLOUD2_SHAPE),0.0);
  //u = 3.0*u*u - 2.0*u*u*u;
  vec2 v = 1.0 - u;

  // rain transition
  vec2 t = vec2(0.2001+0.2*rain, 0.1999+0.2*rain*rain);

  // mix noise gradients
  float n = v.y*(randt(p0,t)*v.x + randt(p0+vec2(1.0,0.0),t)*u.x) +
            u.y*(randt(p0+vec2(0.0,1.0),t)*v.x + randt(p0+vec2(1.0,1.0),t)*u.x);

  float b = bevel(2.0*abs(pos.y-0.5), 0.3);
  return n*b;
}

vec4 renderClouds(vec3 vDir, vec3 vPos, float rain, float time, vec3 fog_col, vec3 sky_col) {
  // local cloud pos
  vec3 pos = vPos;
  pos.y = 0.0;
  pos.xz = NL_CLOUD2_SCALE*(vPos.xz + vec2(1.0,0.5)*(time*NL_CLOUD2_VELOCIY));

  // scaled ray offset
  float height = 7.0*(NL_CLOUD2_THICKNESS + rain*(NL_CLOUD2_RAIN_THICKNESS - NL_CLOUD2_THICKNESS));
  vec3 delta_p;
  delta_p.xz = (NL_CLOUD2_SCALE*height)*vDir.xz;
  delta_p.y = vDir.y;
  delta_p /= float(NL_CLOUD2_STEPS)*(0.02+0.98*abs(vDir.y));

  // alpha, gradient, ray depth temp
  vec3 d = vec3(0.0,1.0,1.0);
  for (int i=0; i<NL_CLOUD2_STEPS; i++) {
    pos += delta_p;
    float m = cloud_df(pos.xyz, rain);
    d.x += m;
    d.y = mix(d.y, pos.y, d.z);
    d.z *= 1.0 - m;
  }
  d.x = d.x/(float(NL_CLOUD2_STEPS)/NL_CLOUD2_DENSITY + d.x);
  d.x = smoothstep(0.2,0.4,d.x);
  if (vPos.y > 0.0) {
  	     d.y = 1.0 - d.y;
  }

  d.y = 1.0-0.8*d.y*d.y;

  vec4 col = vec4(0.0*sky_col, d.x);
  col.rgb += (vec3(0.03,0.05,0.05) + 1.23*fog_col)*d.y;
  col.rgb *= 1.0 - 0.5*rain;

  return col;
}

// aurora is rendered on clouds layer
#ifdef NL_AURORA
vec4 renderAurora(vec3 p, float t, float rain, vec3 sky_col) {
  p.xz *= NL_AURORA_SCALE;
  t *= NL_AURORA_VELOCITY;
  p.xz += 0.05*sin(p.x*4.0 + 20.0*t);
  float d0 = sin(p.x*0.1 + t + sin(p.z*0.2));
  float d1 = sin(p.z*0.1 - t + sin(p.x*0.2));
  float d2 = sin(p.z*0.1 + 1.0*sin(d0 + d1*2.0) + d1*2.0 + d0*1.0);
  d0 *= d0; d1 *= d1; d2 *= d2;
  d2 = d0/(1.0 + d2/NL_AURORA_WIDTH);

  float mask = (1.0-0.8*rain)/(1.0 + 64.0*sky_col.b*sky_col.b);
  return vec4(NL_AURORA*mix(NL_AURORA_COL1,NL_AURORA_COL2,d1),1.0)*d2*mask;
}
#endif

#endif
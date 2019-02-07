#version 300 es
precision highp float;

// The vertex shader used to render the background of the scene

in vec4 vs_Pos;
out vec2 fs_Pos;
const vec4 lightPos = vec4(0, 15, 10, 1);
out vec4 fs_LightVec;

void main() {
  fs_Pos = vs_Pos.xy;
  fs_LightVec = lightPos;
  gl_Position = vs_Pos;
}

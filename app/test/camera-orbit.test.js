import test from 'node:test';
import assert from 'node:assert/strict';
import * as THREE from 'three';
import {mountCameraOrbit} from '../src/camera-orbit.js';

// A DOM stub with the two things the orbit reads -- a height in pixels, so a
// drag can be converted into an angle, and pointer events -- and nothing else.
function element(height=1000){
 const listeners=new Map();
 return {clientHeight:height,clientWidth:1440,style:{},
  addEventListener:(type,fn)=>{(listeners.get(type)??listeners.set(type,[]).get(type)).push(fn);},
  removeEventListener:(type,fn)=>{const list=listeners.get(type)||[];const i=list.indexOf(fn);if(i>=0)list.splice(i,1);},
  setPointerCapture(){},releasePointerCapture(){},
  emit(type,event){for(const fn of listeners.get(type)||[])fn({pointerId:1,button:0,pointerType:'mouse',preventDefault(){},...event});}};
}
function rig({position=[0,0,4],up=[0,1,0],height=1000}={}){
 const camera=new THREE.PerspectiveCamera(35,1.44,.01,1000);
 camera.position.set(...position);camera.up.set(...up);
 const dom=element(height);
 const controls=mountCameraOrbit(camera,dom);
 controls.update();
 return {camera,dom,controls};
}
// Drags are delivered undamped so one update settles them; damping only spreads
// the same total rotation over frames.
function drag(dom,dx,dy,steps=8){
 dom.emit('pointerdown',{clientX:0,clientY:0});
 for(let i=1;i<=steps;i++)dom.emit('pointermove',{clientX:dx*i/steps,clientY:dy*i/steps});
 dom.emit('pointerup',{clientX:dx,clientY:dy});
}
const direction=(camera,controls)=>camera.position.clone().sub(controls.target).normalize();
const degrees=(a,b)=>THREE.MathUtils.radToDeg(a.angleTo(b));

test('equal drags turn the camera equally, whichever way they go and wherever they start',()=>{
 const starts=[
  {position:[0,0,4],up:[0,1,0]},          // coronal
  {position:[4,0,0],up:[0,1,0]},          // sagittal
  {position:[0,4,0],up:[0,0,-1]},         // transverse, exactly on the old pole
  {position:[0,3.99,.28],up:[0,1,0]},     // four degrees off the superior axis
  {position:[1.9,1.6,2.0],up:[0,1,0]},    // oblique
 ];
 for(const start of starts){
  const readings=[];
  for(const [dx,dy] of [[120,0],[-120,0],[0,120],[0,-120],[85,85],[-85,85]]){
   const {camera,dom,controls}=rig(start);
   const before=direction(camera,controls);
   drag(dom,dx,dy);controls.update();
   readings.push(degrees(before,direction(camera,controls)));
  }
  // One screen height is one full turn, so 120 px of a 1000 px viewport is
  // exactly 43.2 degrees along either screen axis, from every start.
  for(const [i,value] of readings.slice(0,4).entries())
   assert.ok(Math.abs(value-43.2)<.02,
    `${JSON.stringify(start.position)} drag ${i} turned ${value.toFixed(3)} deg`);
  // A diagonal of the same 120.2 px length composes two finite rotations, which
  // SO(3) delivers a shade short of their 43.27 degree sum. That shortfall is
  // geometry, not a privileged axis: it is identical for both diagonals.
  assert.ok(Math.abs(readings[4]-readings[5])<.02);
  assert.ok(Math.abs(readings[4]-43.27)<1.3,`diagonal turned ${readings[4].toFixed(3)} deg`);
  const ratio=Math.max(...readings)/Math.min(...readings);
  assert.ok(ratio<1.03,`anisotropy ${ratio.toFixed(4)} at ${JSON.stringify(start.position)}`);
 }
});

test('a sustained vertical drag carries the camera over the head and round again',()=>{
 const {camera,dom,controls}=rig();
 const start=direction(camera,controls);
 let travelled=0,previous=start.clone(),overhead=false,behind=false;
 for(let i=0;i<12;i++){
  drag(dom,0,-100);controls.update();
  const now=direction(camera,controls);
  travelled+=degrees(previous,now);previous=now.clone();
  if(now.y>.9)overhead=true;
  if(now.z<-.9)behind=true;
 }
 // 1200 px of a 1000 px viewport asks for 432 degrees; a turntable delivers 90
 // and stops on the pole.
 assert.ok(Math.abs(travelled-432)<1,`travelled ${travelled.toFixed(2)} deg`);
 assert.ok(overhead,'the camera never reached the view over the head');
 assert.ok(behind,'the camera never carried on past the head to the back');
 // 432 degrees is one full turn plus 72, so it ends 72 degrees from the start.
 assert.ok(Math.abs(degrees(start,direction(camera,controls))-72)<1);
});

test('drags reach orientations a turntable has no rotation for',()=>{
 // A turntable's camera up is always the superior axis projected into the image
 // plane -- that is the whole of its orientation, and it is why it cannot roll.
 // Four non-commuting drags leave the frame off that plane, which is a rotation
 // about the view axis: at the front, about z, the anterior direction.
 const {camera,dom,controls}=rig();
 controls.rollRelaxFactor=0;   // measure the drags, not the relaxation
 for(const [dx,dy] of [[0,-140],[140,0],[0,140],[-140,0]]){drag(dom,dx,dy);controls.update();}
 const view=direction(camera,controls);
 const turntableUp=new THREE.Vector3(0,1,0).addScaledVector(view,-view.y).normalize();
 const lean=Math.min(camera.up.angleTo(turntableUp),camera.up.angleTo(turntableUp.clone().negate()));
 assert.ok(THREE.MathUtils.radToDeg(lean)>10,
  `expected a roll no turntable can make, got ${THREE.MathUtils.radToDeg(lean).toFixed(2)} deg`);
});

test('an idle view relaxes its lean back toward the superior axis without flipping',()=>{
 const {camera,dom,controls}=rig();
 // Roll the front view 30 degrees, the lean a non-commuting pair of drags leaves.
 camera.up.applyAxisAngle(new THREE.Vector3(0,0,1),THREE.MathUtils.degToRad(30));
 controls.update();
 for(let i=0;i<400;i++)controls.update();
 assert.ok(camera.up.angleTo(new THREE.Vector3(0,1,0))<1e-3,
  `lean left at ${THREE.MathUtils.radToDeg(camera.up.angleTo(new THREE.Vector3(0,1,0))).toFixed(3)} deg`);
 assert.ok(degrees(direction(camera,controls),new THREE.Vector3(0,0,1))<1e-6,
  'relaxing a lean must not move the camera');

 // A view carried over the head is inverted on purpose: relaxing keeps the
 // superior axis vertical in the image but never turns the body back up.
 const back=rig({position:[0,0,-4],up:[0,-1,0]});
 back.camera.up.applyAxisAngle(new THREE.Vector3(0,0,1),THREE.MathUtils.degToRad(20));
 back.controls.update();
 for(let i=0;i<400;i++)back.controls.update();
 assert.ok(back.camera.up.dot(new THREE.Vector3(0,1,0))<0,'the inverted view was flipped upright');
 assert.ok(back.camera.up.angleTo(new THREE.Vector3(0,-1,0))<1e-3,'the inverted view kept its lean');
});

test('a view on the superior axis is stable and still fully steerable',()=>{
 // The transverse snap sits exactly where a turntable degenerates. Nothing here
 // relaxes it, nothing produces a NaN, and a drag still turns it.
 const {camera,dom,controls}=rig({position:[0,4,0],up:[0,0,-1]});
 for(let i=0;i<200;i++)controls.update();
 assert.ok(camera.position.toArray().every(Number.isFinite));
 assert.ok(camera.up.angleTo(new THREE.Vector3(0,0,-1))<1e-6,'the top view rolled on its own');
 const before=direction(camera,controls);
 drag(dom,120,0);controls.update();
 assert.ok(Math.abs(degrees(before,direction(camera,controls))-43.2)<.05);
});

test('damping delivers the whole drag and nothing more',()=>{
 const {camera,dom,controls}=rig();
 controls.enableDamping=true;
 const before=direction(camera,controls);
 drag(dom,120,0);
 for(let i=0;i<600;i++)controls.update();
 assert.ok(Math.abs(degrees(before,direction(camera,controls))-43.2)<.02);
});

test('disabling the controls hands the pointer to the scene gizmo',()=>{
 const {camera,dom,controls}=rig();
 controls.enabled=false;
 const before=camera.position.clone();
 drag(dom,200,140);controls.update();
 assert.ok(camera.position.distanceTo(before)<1e-9);
});

test('the target, the wheel and a right drag behave as the rest of the app expects',()=>{
 const {camera,dom,controls}=rig();
 controls.target.set(0,.2,0);
 camera.position.set(0,.2,4);
 controls.update();
 const radius=camera.position.distanceTo(controls.target);
 dom.emit('wheel',{deltaY:-120});
 controls.update();
 assert.ok(camera.position.distanceTo(controls.target)<radius,'scrolling up must move closer');
 assert.ok(direction(camera,controls).angleTo(new THREE.Vector3(0,0,1))<1e-6,'a dolly must not turn the camera');
 const target=controls.target.clone();
 dom.emit('pointerdown',{button:2,clientX:0,clientY:0});
 dom.emit('pointermove',{button:2,clientX:60,clientY:0});
 dom.emit('pointerup',{button:2,clientX:60,clientY:0});
 controls.update();
 assert.ok(controls.target.distanceTo(target)>1e-3,'a right drag must pan the target');
 assert.ok(Math.abs(controls.target.y-target.y)<1e-9,'a horizontal pan must not move the target vertically');
});

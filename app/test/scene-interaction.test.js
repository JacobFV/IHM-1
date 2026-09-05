import test from 'node:test';
import assert from 'node:assert/strict';
import {cursorSpring,rotateOffset,advanceScene} from '../src/scene-forces.js';

test('cursor spring uses displacement, preserves direction and bounds force',()=>{
  assert.deepEqual(cursorSpring([0,0,0],[.1,0,0]),[8,0,0]);
  assert.deepEqual(cursorSpring([1,2,3],[1,2,3]),[0,0,0]);
  const f=cursorSpring([0,0,0],[3,4,0]);assert.ok(Math.abs(Math.hypot(...f)-100)<1e-10);assert.equal(f[0]/f[1],3/4);
  assert.throws(()=>cursorSpring([0,0,0],[NaN,0,0]));
});

test('grab offsets follow sphere material rotation between commands',()=>{
  const rotation=[[0,-1,0],[1,0,0],[0,0,1]],point=[.05,0,0];
  assert.deepEqual(rotateOffset(point,rotation),[0,.05,0]);
  assert.deepEqual(rotateOffset(rotateOffset(point,rotation),rotation,true),point);
});

test('uncertain advance reads current state without repeating force',async()=>{
  const calls=[],committed={sequence:1,time_s:.02};
  const result=await advanceScene(async(path,data)=>{
    calls.push({path,data});if(data)throw Error('Response lost');return committed;
  },'/session',{sequence:0,seconds:.02});
  assert.equal(result.frame,committed);assert.equal(result.recovered,true);
  assert.deepEqual(calls.map(c=>c.path),['/session/step','/session']);
  await assert.rejects(advanceScene(async()=>{throw Error('Offline');},'/session',{}),/Offline/);
});

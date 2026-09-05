import test from 'node:test';
import assert from 'node:assert/strict';
import {ElasticHairSystem,createHairController} from '../src/hair_dynamics.js';
const fixture=()=>({centerlines_m:[0,0,0,.005,0,0,.01,0,0,.015,0,0,.02,0,0,.025,0,0,.03,0,0],strand_offsets:[0,7],radius_m:[40e-6],tensile_modulus_pa:7.11e9,bending_modulus_pa:5.7e9,density_kg_m3:1312});
test('physical strand mass and exact tube radius',()=>{
 const s=new ElasticHairSystem(fixture());assert.ok(Math.abs(s.totalMassKg()-1312*Math.PI*(40e-6)**2*.03)<1e-20);
 const tube=s.writeTubePositions();assert.ok(Math.abs(Math.hypot(tube[0],tube[1],tube[2])-40e-6)<1e-10);
});
test('gravity bends a rooted strand rather than transporting a rigid shaft',()=>{
 const s=new ElasticHairSystem(fixture(),{maxSubstep:1/1000});
 for(let k=0;k<50;k++)s.step(.001);assert.deepEqual(Array.from(s.positions.slice(0,6)),[0,0,0,.005,0,0]);
 assert.ok(s.positions[19]<-1e-6);assert.ok(s.maxStrain()<.01);assert.ok(Number.isFinite(s.energy()));
});
test('root translation is exact and root work does not masquerade as gravity',()=>{
 const s=new ElasticHairSystem(fixture());s.step(.001,{roots_m:[.1,.2,.3],root_tangents:[0,1,0],gravity_m_s2:[0,0,0]});
 assert.deepEqual(Array.from(s.positions.slice(0,3)),[.1,.2,.3]);assert.ok(Math.abs(s.positions[4]-.205)<1e-15);
});
test('stationary zero-force state and checkpoint replay are deterministic',()=>{
 const s=new ElasticHairSystem(fixture());const start=s.positions.slice();s.step(.01,{gravity_m_s2:[0,0,0]});assert.deepEqual(s.positions,start);
 s.step(.005);const checkpoint=s.checkpoint();s.step(.005);const end=s.positions.slice();s.restore(checkpoint);s.step(.005);assert.deepEqual(s.positions,end);
});
test('invalid dt is isolated before changing state',()=>{
 const s=new ElasticHairSystem(fixture());const before=s.checkpoint();assert.throws(()=>s.step(-1));assert.throws(()=>s.step(.1,{gravity_m_s2:[NaN,0,0]}));assert.deepEqual(s.checkpoint(),before);
});
test('fixed-root gravity integration does not create mechanical energy',()=>{
 const s=new ElasticHairSystem(fixture(),{maxSubstep:1/240});
 for(let k=0;k<24;k++){s.step(1/240);assert.ok(s.energy()<=1e-15,'Unforced fixed-root mechanical energy increased');}
});
test('timestep refinement bounds the short-hair transient',()=>{
 const tips=[];for(const dt of[1/240,1/480]){const s=new ElasticHairSystem(fixture(),{maxSubstep:dt});for(let k=0;k<Math.round(.1/dt);k++)s.step(dt);tips.push(s.positions[19]);}
 assert.ok(Math.abs(tips[0]-tips[1])<15e-6);
});
test('external tip force causes a distinct elastic response',()=>{
 const a=new ElasticHairSystem(fixture()),b=new ElasticHairSystem(fixture()),force=new Float64Array(21);force[20]=1e-6;
 for(let k=0;k<24;k++){a.step(1/240);b.step(1/240,{external_forces_n:force});}
 assert.ok(b.positions[20]>1e-5);assert.equal(a.positions[20],0);
});

const controllerFixture=()=>({strands:fixture(),attachment:{skin_entity_id:'skin',barycentric:[[1,0,0]],reference_triangles_m:[0,0,0,0,1,0,0,0,1]}});
test('controller limits updates to 20 Hz, pauses and explicitly resets missing time',()=>{
 const c=createHairController(controllerFixture(),{maxUpdateHz:20});assert.equal(c.update({simulationTime:0}).diagnostics.reset_reason,'initial');
 assert.equal(c.update({simulationTime:.01}).diagnostics.updated,false);assert.equal(c.system.time,0);
 assert.equal(c.update({simulationTime:.05}).diagnostics.updated,true);assert.equal(c.system.time,.05);
 assert.equal(c.update({simulationTime:.1,visible:false}).diagnostics.paused,true);assert.equal(c.system.time,.05);
 assert.equal(c.update({simulationTime:.2}).diagnostics.reset_reason,'resume_without_catchup');
 assert.equal(c.update({simulationTime:5}).diagnostics.reset_reason,'time_gap_not_simulated');
 assert.equal(c.update({simulationTime:0}).diagnostics.reset_reason,'rewind');
});
test('render fibers inherit computed guide deformation while roots remain attached',()=>{
 const g=controllerFixture();g.render_strands=fixture();g.render_attachment=g.attachment;g.guide_interpolation=[[[0,1]]];const c=createHairController(g,{maxUpdateHz:20});
 const start=c.update({simulationTime:0}).positions.slice();const moved=c.update({simulationTime:.05});assert.equal(moved.diagnostics.simulated_guides,1);assert.equal(moved.diagnostics.rendered_fibers,1);assert.notDeepEqual(moved.positions,start);assert.deepEqual(moved.positions.slice(0,12),start.slice(0,12));
});

import * as THREE from 'three';
import { TransformControls } from 'three/addons/controls/TransformControls.js';
import './scene-interaction.css';
import {cursorSpring,advanceScene} from './scene-forces.js';
import {bodyEndpoint,bodyEnvironment,bodyCommand,frameScope,materialOffset,materialPoint,createBodyOwner,closeBodyOwner,scheduleBodyIntakes} from './embodied-live.js';
import {mountTemporalSpectrumMonitor} from './temporal-spectrum-monitor.js';
import {mountSkinVoltageMonitor} from './skin-voltage-monitor.js';
import {mountIntakeMassMonitor} from './intake-mass-monitor.js';
import {mountIntakeMonitor} from './intake-monitor.js';
import {mountEmbodiedPanels} from './embodied-panels.js';

export function mountSceneInteraction({scene,camera,renderer,controls,group,getObjects,getSelected,onSelect,onFrame,onPauseReplay}) {
  const mount=document.getElementById('scene-controls'),monitor=document.getElementById('scene-monitor');
  if(!mount||!monitor)throw Error('Scene interaction mounts missing');
  mount.innerHTML=`<label class="field-label" for="scene-environment">Environment</label>
    <select id="scene-environment"><option value="bed">Bed · supine contact</option><option value="floor">Floor · upright contact</option><option value="studio">Free space · no gravity</option></select>
    <div class="scene-modes" role="group" aria-label="Cursor mode"><button data-scene-mode="select" aria-pressed="true">Select</button><button data-scene-mode="gimbal" aria-pressed="false">Gimbal</button><button data-scene-mode="force" aria-pressed="false">Force</button></div>
    <div class="scene-actions"><button id="scene-play">Start Body</button><button id="scene-reset">Reset</button></div>
    <p id="scene-status-note" class="muted" role="status">Select a part to inspect. Start Body initializes the unified native body.</p><details class="scene-advanced"><summary>Advanced execution</summary><label class="field-label" for="scene-owner">Execution owner</label><button id="scene-reconnect" type="button">Reconnect existing body</button><select id="scene-owner"><option value="embodied">Unified body · native + neural</option><option value="reduced">Reduced mechanics experiment</option></select><p class="muted">The reduced experiment has no coupled native physiology. No automatic fallback occurs.</p></details>`;
  monitor.innerHTML=`<div class="scene-values"><span>Time <b id="scene-clock">0.00 s</b></span><span>Applied force <b id="scene-force-value">0 N</b></span></div><p id="scene-target" class="muted">No force target</p><details class="evidence-fold"><summary>Mechanics and supports</summary><p id="scene-scope">The unified body advances articulation, muscle sensors, pinned neural dynamics and physiology together. Calibration and full tissue/garment integration remain incomplete.</p></details>`;
  const $=id=>document.getElementById(id);
  const environmentGroup=new THREE.Group();environmentGroup.name='Interactive environment';group.add(environmentGroup);
  const objectMeshes=new Map(),ray=new THREE.Raycaster(),plane=new THREE.Plane(),cursor=new THREE.Vector3();
  const gizmoTarget=new THREE.Object3D();group.add(gizmoTarget);
  const gizmo=new TransformControls(camera,renderer.domElement);gizmo.setMode('translate');gizmo.setSize(.65);
  scene.add(gizmo.getHelper());
  const arrow=new THREE.ArrowHelper(new THREE.Vector3(0,1,0),new THREE.Vector3(),0,0xe9b979,.02,.009);group.add(arrow);arrow.visible=false;
  let mode='select',session=null,state=null,running=false,pending=false,creating=null,drag=null,gizmoDragging=false;
  let lastStep=0,lastPoll=0,disposed=false,environment='bed',kind='embodied',initializing=false,faulted=false,selectedId=null,lastError='',resetting=false,resetTask=null;
  const panels=mountEmbodiedPanels();
  const skinVoltage=mountSkinVoltageMonitor(document.getElementById('skin-voltage-monitor'));
  const intakeMass=mountIntakeMassMonitor(document.getElementById('intake-mass-monitor'));
  const spectrum=mountTemporalSpectrumMonitor(document.getElementById('temporal-spectrum-monitor'));
  let intakeTask=null,intakeOwner=null,intakeSignature='',intakeConnected=false;
  let intake=mountIntakeMonitor(document.getElementById('intake-monitor'),{submit:scheduleIntakes});
  function syncIntake(){const connected=!disposed&&!resetting&&!initializing&&!faulted&&kind==='embodied'&&!!session&&!!state;if(connected!==intakeConnected){intakeConnected=connected;intake.setConnected(connected);}}
  function clearIntake(){intake.dispose();intake=mountIntakeMonitor(document.getElementById('intake-monitor'),{submit:scheduleIntakes});intakeOwner=null;intakeSignature='';intakeConnected=false;}
  async function scheduleIntakes(events) {
    if(intakeTask)throw Error('An intake request is already pending');
    if(disposed||resetting||initializing||faulted||kind!=='embodied'||!session||!state)throw Error('Start a unified body before scheduling intake');
    const owner=session;
    intakeTask=Promise.resolve().then(async()=>{
      while(pending)await new Promise(resolve=>setTimeout(resolve,10));
      if(disposed||resetting||owner!==session)throw Error('Body owner changed before intake could be scheduled');
      pending=true;
      try{
        const result=await scheduleBodyIntakes(request,endpoint()+'/'+owner,state.sequence,events);
        if(disposed)return;
        accept(result.frame);
        if(result.recovered){
          running=false;playText('Resume '+label());
          if(result.error.definitelyRejected)throw result.error;
          const known=new Set((result.frame.intake_schedule?.events||[]).map(e=>e.event_id));
          if(!events.every(e=>known.has(e.event_id)))throw result.error;
          status('Intake request recovered from current body state. Body paused; review the schedule before resuming.');
        }
        return result.frame.intake_schedule;
      }catch(error){running=false;playText(faulted?'Reset required':'Resume '+label());status(error.message);throw error;}
      finally{pending=false;}
    });
    try{return await intakeTask;}finally{intakeTask=null;syncIntake();}
  }
  const endpoint=()=>bodyEndpoint(kind);
  const playText=(text)=>{if(!disposed)$('scene-play').textContent=text;};
  const label=()=>kind==='embodied'?'Body':'reduced mechanics';
  const status=text=>{if(!disposed)$('scene-status-note').textContent=text;};

  async function request(path,data) {
    const response=await fetch(path,data===undefined?{cache:'no-store'}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
    const payload=await response.json();if(!response.ok){const error=Error(payload.error||'Scene request failed');error.httpStatus=response.status;throw error;}return payload;
  }
  function allObjects() {
    const raw=getObjects();return [...(raw instanceof Map?raw.values():raw||[]),...objectMeshes.values()].filter(x=>x.visible);
  }
  function pick(event) {
    const r=renderer.domElement.getBoundingClientRect();ray.setFromCamera(new THREE.Vector2((event.clientX-r.left)/r.width*2-1,-(event.clientY-r.top)/r.height*2+1),camera);
    const hit=ray.intersectObjects(allObjects(),true)[0];if(!hit)return null;
    let mesh=hit.object;while(!mesh.userData.structure&&!mesh.userData.sceneObject&&mesh.parent)mesh=mesh.parent;
    const ident=mesh.userData.sceneObject||mesh.userData.structure?.id;if(!ident)return null;
    return {id:ident,point:group.worldToLocal(hit.point.clone()),mesh,world:hit.point.clone()};
  }
  function centroid(id,fallback) {
    const object=state?.objects?.find(o=>o.id===id),body=state?.entities?.[id];
    if(object)return new THREE.Vector3(...object.position_m);
    if(body)return new THREE.Vector3(...body.centroid_m);
    if(fallback?.geometry){fallback.geometry.computeBoundingBox();const p=fallback.geometry.boundingBox.getCenter(new THREE.Vector3());return group.worldToLocal(fallback.localToWorld(p));}
    return new THREE.Vector3();
  }
  function selected(hit) {
    selectedId=hit.id;if(hit.mesh.userData.structure)onSelect(hit.mesh.userData.structure);
    $('scene-target').textContent=hit.mesh.userData.structure?.name||'Free sphere · 0.4 kg';
    gizmoTarget.position.copy(centroid(hit.id,hit.mesh));gizmoTarget.updateMatrixWorld();
    if(mode==='gimbal')gizmo.attach(gizmoTarget);else gizmo.detach();
  }
  function environmentMeshes() {
    while(environmentGroup.children.length){const o=environmentGroup.children[0];environmentGroup.remove(o);o.geometry?.dispose();o.material?.dispose();}
    objectMeshes.clear();
    if(kind==='reduced'&&environment!=='studio'){
      const bed=environment==='bed';
      const mesh=new THREE.Mesh(new THREE.BoxGeometry(bed?1.05:3,bed?2.1:.06,bed?.08:3),new THREE.MeshStandardMaterial({color:bed?0x718680:0x344646,roughness:.95}));
      if(bed)mesh.position.set(0,0,-.28);else mesh.position.set(0,-.99,0);
      environmentGroup.add(mesh);
    }
  }
  function accept(frame) {
    if(!frame?.entities||!Number.isFinite(frame.time_s)||!Number.isInteger(frame.sequence)||(kind==='embodied'&&frame.schema!=='ihm.embodied-frame.v1')){faulted=true;syncIntake();throw Error(frame?.error||'Body owner has no valid live frame; Reset required.');}
    panels.update(session,frame);state=frame;intakeMass.update(session,frame);skinVoltage.update(session,frame);
    if(kind==='embodied')spectrum.update(session,frame);else spectrum.clear();
    if(kind==='embodied'){
      if(intakeOwner&&intakeOwner!==session)clearIntake();intakeOwner=session;
      const schedule=frame.intake_schedule||{events:[]},signature=JSON.stringify(schedule);
      if(signature!==intakeSignature){intake.update(schedule);intakeSignature=signature;}
    }
    syncIntake();onFrame(frame);
    const currentEnvironment={free:'studio',supine:'bed',upright:'floor'}[frame.mechanics?.body_environment?.kind];
    if(kind==='embodied'&&currentEnvironment){environment=currentEnvironment;$('scene-environment').value=environment;}
    $('scene-scope').textContent=frameScope(frame);
    for(const item of frame.objects||[]){
      let mesh=objectMeshes.get(item.id);
      if(!mesh){mesh=new THREE.Mesh(new THREE.SphereGeometry(item.radius_m,24,16),new THREE.MeshStandardMaterial({color:0xd6a16e,roughness:.6}));mesh.userData.sceneObject=item.id;environmentGroup.add(mesh);objectMeshes.set(item.id,mesh);}
      mesh.position.fromArray(item.position_m);const r=item.rotation_matrix;mesh.setRotationFromMatrix(new THREE.Matrix4().set(r[0][0],r[0][1],r[0][2],0,r[1][0],r[1][1],r[1][2],0,r[2][0],r[2][1],r[2][2],0,0,0,0,1));
    }
    $('scene-clock').textContent=frame.time_s.toFixed(2)+' s';
    if(selectedId&&!gizmoDragging)gizmoTarget.position.copy(centroid(selectedId));
  }
  function initialized(frame) {
    initializing=false;onPauseReplay();accept(frame);
    if(drag?.initialAnchor){
      if(drag.gimbal)drag.target.sub(drag.initialAnchor).add(centroid(drag.id));
      else drag.offset.fromArray(materialOffset(drag.initialAnchor.toArray(),materialEntity(drag.id)));
      delete drag.initialAnchor;
    }
    status(kind==='embodied'?'Unified body ready · live inputs act on the next body tick.':'Reduced mechanics experiment ready · no native physiology feedback.');
    playText(running?'Pause '+label():'Resume '+label());
  }
  function pendingStatus(frame) {
    if(frame?.schema==='ihm.embodied-frame.v1'){initialized(frame);return;}
    if(frame?.error||frame?.closed||frame?.status==='error')throw Error(frame.error||'Body closed before initialization completed.');
    if(!['initializing','ready','closing'].includes(frame?.status))throw Error('Unexpected body startup state');
    initializing=true;syncIntake();status('Body initializing · waiting for the native resource slot. Reset requests cleanup.');
    panels.status('Body initializing · no live physiological frame yet.');
  }
  async function start() {
    if(disposed||resetting)return;
    if(faulted){status('Reset is required before this body can restart.');return;}
    if(creating)return creating;
    if(!session){
      creating=(async()=>{
        onPauseReplay();running=true;playText('Initializing '+label()+'…');status('Initializing '+label()+'…');
        let frame;
        if(kind==='embodied'){
          frame=await createBodyOwner(request,endpoint(),{environment:bodyEnvironment(environment,kind)});
        }else frame=await request(endpoint(),{environment:bodyEnvironment(environment,kind)});
        session=frame.id;
        if(!session)throw Error('Body startup returned no session identity');
        if(disposed){await request(endpoint()+'/'+session+'/close',{});return;}
        if(kind==='embodied')pendingStatus(frame);else initialized(frame);
      })();
      try{await creating;}catch(error){running=false;if(session)faulted=true;playText(session?'Reset required':'Start '+label());throw error;}finally{creating=null;}
    }else if(!initializing){running=true;playText('Pause '+label());panels.status('Body resumed · awaiting the next accepted native frame.');}
  }
  function reset() {
    if(resetTask)return resetTask;
    resetting=true;syncIntake();running=false;drag=null;gizmoDragging=false;arrow.visible=false;gizmo.detach();controls.enabled=true;
    resetTask=Promise.resolve().then(async()=>{
    try{
      if(creating){try{await creating;}catch(error){if(!session)throw error;}}
      while(pending||intakeTask)await new Promise(resolve=>setTimeout(resolve,10));
      running=false;
      if(session){
        if(kind==='embodied')await closeBodyOwner(request,endpoint()+'/'+session,async()=>{
          status('Closing body · waiting for native cleanup before changing views.');panels.status('Closing body · last accepted state remains shown until cleanup completes.');
          await new Promise(resolve=>setTimeout(resolve,500));
        });
        else if((await request(endpoint()+'/'+session+'/close',{})).closed!==true)throw Error('Reduced scene did not confirm closure');
      }
      if(disposed)return;
      session=null;state=null;initializing=false;faulted=false;selectedId=null;environment=$('scene-environment').value;kind=$('scene-owner').value;environmentMeshes();panels.clear();spectrum.clear();intakeMass.clear();skinVoltage.clear();clearIntake();syncIntake();
      $('scene-play').textContent='Start '+label();$('scene-clock').textContent='0.00 s';$('scene-force-value').textContent='0 N';
      onFrame(null);status('Body reset. Start '+label()+' to advance.');
    }catch(e){status(e.message);throw e;}finally{resetting=false;resetTask=null;syncIntake();}
    });
    return resetTask;
  }
  $('scene-reconnect').onclick=async()=>{
    if(disposed||resetting||pending||intakeTask||creating)return;
    if(kind!=='embodied'){status('Select Unified body before reconnecting.');return;}
    running=false;pending=true;
    try{
      const list=await request(endpoint());const active=list.sessions.filter(s=>!s.closed);
      if(active.length!==1)throw Error(active.length?'More than one body owner exists; no automatic choice made.':'No existing live body is available.');
      session=active[0].id;const frame=await request(endpoint()+'/'+session);if(disposed||resetting)return;faulted=false;onPauseReplay();pendingStatus(frame);
      status('Existing body reconnected and paused. Resume only from this controlling view.');
    }catch(error){status(error.message);}finally{pending=false;}
  };
  $('scene-owner').onchange=()=>reset().catch(e=>status(e.message));
  $('scene-environment').onchange=()=>reset().catch(e=>status(e.message));$('scene-reset').onclick=()=>reset().catch(e=>status(e.message));
  $('scene-play').onclick=()=>{if(running){running=false;$('scene-play').textContent='Resume '+label();panels.status('Paused · last accepted body state. Inputs take effect on resume.');}else start().catch(e=>status(e.message));};
  mount.querySelectorAll('[data-scene-mode]').forEach(button=>button.onclick=()=>{
    mode=button.dataset.sceneMode;mount.querySelectorAll('[data-scene-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b===button)));
    drag=null;arrow.visible=false;gizmo.detach();controls.enabled=true;
    if(mode==='gimbal'&&selectedId)gizmo.attach(gizmoTarget);
    status(mode==='select'?'Click a structure to inspect.':mode==='gimbal'?'Select a part; drag an axis to apply a cursor spring.':'Drag a body part or free object to apply force.');
  });
  gizmo.addEventListener('dragging-changed',event=>{
    gizmoDragging=event.value;controls.enabled=!event.value;
    if(event.value&&selectedId){drag={id:selectedId,offset:new THREE.Vector3(),target:gizmoTarget.position.clone(),initialAnchor:session&&!initializing?null:gizmoTarget.position.clone(),gimbal:true};start().catch(e=>{status(e.message);gizmoDragging=false;drag=null;});}
    else {drag=null;arrow.visible=false;}
  });
  gizmo.addEventListener('objectChange',()=>{if(gizmoDragging&&drag)drag.target.copy(gizmoTarget.position);});
  function down(event) {
    if(mode==='select'||event.button!==0||gizmoDragging||resetting)return;
    // Axis handles own their pointer events; only pick new anatomy away from an axis.
    if(mode==='gimbal'&&gizmo.axis)return;
    const hit=pick(event);if(!hit)return;selected(hit);
    if(mode==='gimbal')return;
    let offset;try{offset=session&&!initializing?materialOffset(hit.point.toArray(),materialEntity(hit.id)):[0,0,0];}catch(error){status(error.message);return;}
    event.preventDefault();event.stopImmediatePropagation();controls.enabled=false;
    renderer.domElement.setPointerCapture(event.pointerId);
    plane.setFromNormalAndCoplanarPoint(camera.getWorldDirection(new THREE.Vector3()),hit.world);
    drag={id:hit.id,offset:new THREE.Vector3(...offset),target:hit.point.clone(),initialAnchor:session&&!initializing?null:hit.point.clone(),pointer:event.pointerId};
    start().catch(e=>{status(e.message);drag=null;controls.enabled=true;});
  }
  function move(event) {
    if(!drag||mode!=='force'||drag.pointer!==event.pointerId)return;
    event.preventDefault();event.stopImmediatePropagation();
    const r=renderer.domElement.getBoundingClientRect();ray.setFromCamera(new THREE.Vector2((event.clientX-r.left)/r.width*2-1,-(event.clientY-r.top)/r.height*2+1),camera);
    if(ray.ray.intersectPlane(plane,cursor))drag.target.copy(group.worldToLocal(cursor.clone()));
  }
  function up(event) {
    if(mode!=='force'||!drag)return;
    event.preventDefault();event.stopImmediatePropagation();drag=null;controls.enabled=true;arrow.visible=false;$('scene-force-value').textContent='0 N';
    if(renderer.domElement.hasPointerCapture(event.pointerId))renderer.domElement.releasePointerCapture(event.pointerId);
  }
  renderer.domElement.addEventListener('pointerdown',down,true);renderer.domElement.addEventListener('pointermove',move,true);
  renderer.domElement.addEventListener('pointerup',up,true);renderer.domElement.addEventListener('pointercancel',up,true);
  function materialEntity(id) {
    const entity=state?.entities?.[id];if(entity)return entity;
    const object=state?.objects?.find(o=>o.id===id);
    if(object)return {...object,centroid_m:object.position_m};
    throw Error('Selected structure has no live mechanical owner');
  }
  function forceCommand() {
    if(!drag)return [];
    const point=new THREE.Vector3(...materialPoint(drag.offset.toArray(),materialEntity(drag.id))),force=cursorSpring(point.toArray(),drag.target.toArray());
    const magnitude=Math.hypot(...force);$('scene-force-value').textContent=magnitude.toFixed(2)+' N';
    arrow.position.copy(point);arrow.visible=magnitude>1e-5;
    if(arrow.visible){arrow.setDirection(new THREE.Vector3(...force).normalize());arrow.setLength(Math.min(.3,magnitude*.015),.018,.009);}
    return [{id:drag.id,force_n:force,point_m:point.toArray()}];
  }
  async function update(now) {
    if(disposed||resetting||pending||intakeTask||!session)return;
    if(initializing){
      if(now-lastPoll<500)return;lastPoll=now;pending=true;
      try{const frame=await request(endpoint()+'/'+session);if(!disposed)pendingStatus(frame);}
      catch(error){running=false;initializing=false;faulted=true;playText('Reset required');status(error.message);if(!disposed)panels.status('Body startup stopped: '+error.message);}
      finally{pending=false;}
      return;
    }
    if(!running||now-lastStep<20)return;
    pending=true;lastStep=now;
    try {
      const result=await advanceScene(request,endpoint()+'/'+session,bodyCommand(kind,state.sequence,forceCommand(),panels.inputs,state));
      if(!disposed)accept(result.frame);
      if(result.recovered)throw Error(`${result.error.message}. Recorded state restored; resume when ready.`);
    }
    catch(e){running=false;drag=null;arrow.visible=false;controls.enabled=true;if(!disposed){panels.status('Stopped · last accepted body frame. '+e.message);$('scene-play').textContent=faulted?'Reset required':'Resume '+label();if(e.message!==lastError)status(e.message);}lastError=e.message;}
    finally{pending=false;}
  }
  function unload(){if(session)navigator.sendBeacon(endpoint()+'/'+session+'/close',new Blob(['{}'],{type:'application/json'}));}
  window.addEventListener('pagehide',unload);environmentMeshes();
  return {update,reset,scheduleIntakes,createControls:()=>mount,get state(){return state;},get active(){return !!session||!!creating||pending;},
    dispose(){disposed=true;skinVoltage.dispose();intakeMass.dispose();spectrum.dispose();intake.dispose();running=false;unload();window.removeEventListener('pagehide',unload);gizmo.dispose();scene.remove(gizmo.getHelper());group.remove(environmentGroup,gizmoTarget,arrow);
      renderer.domElement.removeEventListener('pointerdown',down,true);renderer.domElement.removeEventListener('pointermove',move,true);renderer.domElement.removeEventListener('pointerup',up,true);renderer.domElement.removeEventListener('pointercancel',up,true);}};
}

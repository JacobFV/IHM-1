import * as THREE from 'three';
import { TransformControls } from 'three/addons/controls/TransformControls.js';
import {cursorSpring,advanceScene} from './scene-forces.js';
import {bodyEndpoint,bodyEnvironment,bodyCommand,frameScope,materialOffset,materialPoint,createBodyOwner,closeBodyOwner,scheduleBodyIntakes} from './embodied-live.js';
import {mountTemporalSpectrumMonitor} from './temporal-spectrum-monitor.js';
import {mountSkinVoltageMonitor} from './skin-voltage-monitor.js';
import {mountIntakeMassMonitor} from './intake-mass-monitor.js';
import {mountIntakeMonitor} from './intake-monitor.js';
import {mountEmbodiedPanels} from './embodied-panels.js';

// One native body. There is no owner picker, no candidate picker and no
// reconnect flow: the workbench resolves the active source itself.
export function mountSceneInteraction({scene,camera,renderer,controls,group,getObjects,onSelect,onFrame,onPauseReplay,mount,monitor,onStatus}) {
  mount.innerHTML=`<div class="scene-modes" role="group" aria-label="Cursor mode"><button data-scene-mode="select" aria-pressed="true">Select</button><button data-scene-mode="gimbal" aria-pressed="false">Gimbal</button><button data-scene-mode="force" aria-pressed="false">Force</button></div><div class="scene-actions"><button id="scene-reset" type="button">Reset</button></div>`;
  monitor.innerHTML=`<div class="scene-values"><span>Time <b id="scene-clock">0.00 s</b></span><span>Applied force <b id="scene-force-value">0 N</b></span></div><p id="scene-target" class="note">No force target</p><p id="scene-scope" class="note"></p>`;
  const $=id=>document.getElementById(id);
  const objectMeshes=new Map(),ray=new THREE.Raycaster(),plane=new THREE.Plane(),cursor=new THREE.Vector3();
  const environmentGroup=new THREE.Group();environmentGroup.name='Interactive environment';group.add(environmentGroup);
  const gizmoTarget=new THREE.Object3D();group.add(gizmoTarget);
  const gizmo=new TransformControls(camera,renderer.domElement);gizmo.setMode('translate');gizmo.setSize(.65);
  scene.add(gizmo.getHelper());
  const arrow=new THREE.ArrowHelper(new THREE.Vector3(0,1,0),new THREE.Vector3(),0,0xe9b979,.02,.009);group.add(arrow);arrow.visible=false;
  let mode='select',session=null,state=null,running=false,pending=false,creating=null,drag=null,gizmoDragging=false;
  let lastStep=0,lastPoll=0,disposed=false,environment='bed',initializing=false,faulted=false,selectedId=null,lastError='',resetting=false,resetTask=null;
  let options={regional_skin:false,intake_mass:false};
  let environmentSelection={};
  const panels=mountEmbodiedPanels();
  const skinVoltage=mountSkinVoltageMonitor($('skin-voltage-monitor'));
  const intakeMass=mountIntakeMassMonitor($('intake-mass-monitor'));
  const spectrum=mountTemporalSpectrumMonitor($('temporal-spectrum-monitor'));
  let intakeTask=null,intakeOwner=null,intakeSignature='',intakeConnected=false;
  let intake=mountIntakeMonitor($('intake-monitor'),{submit:scheduleIntakes});
  const status=text=>{if(!disposed)onStatus(text);};
  function syncIntake(){const connected=!disposed&&!resetting&&!initializing&&!faulted&&!!session&&!!state;if(connected!==intakeConnected){intakeConnected=connected;intake.setConnected(connected);}}
  function clearIntake(){intake.dispose();intake=mountIntakeMonitor($('intake-monitor'),{submit:scheduleIntakes});intakeOwner=null;intakeSignature='';intakeConnected=false;}
  async function scheduleIntakes(events) {
    if(intakeTask)throw Error('An intake request is already pending');
    if(disposed||resetting||initializing||faulted||!session||!state)throw Error('Start the body before scheduling intake');
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
          running=false;
          if(result.error.definitelyRejected)throw result.error;
          const known=new Set((result.frame.intake_schedule?.events||[]).map(e=>e.event_id));
          if(!events.every(e=>known.has(e.event_id)))throw result.error;
          status('Intake request recovered from current body state. Body paused; review the schedule before resuming.');
        }
        return result.frame.intake_schedule;
      }catch(error){running=false;status(error.message);throw error;}
      finally{pending=false;}
    });
    try{return await intakeTask;}finally{intakeTask=null;syncIntake();}
  }
  const endpoint=()=>bodyEndpoint('embodied');

  async function request(path,data) {
    const response=await fetch(path,data===undefined?{cache:'no-store'}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
    const payload=await response.json();if(!response.ok){const error=Error(payload.error||'Body request failed');error.httpStatus=response.status;throw error;}return payload;
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
  function accept(frame) {
    if(!frame?.entities||!Number.isFinite(frame.time_s)||!Number.isInteger(frame.sequence)||frame.schema!=='ihm.embodied-frame.v1'){faulted=true;syncIntake();throw Error(frame?.error||'Body owner has no valid live frame; Reset required.');}
    panels.update(session,frame);state=frame;intakeMass.update(session,frame);skinVoltage.update(session,frame);spectrum.update(session,frame);
    if(intakeOwner&&intakeOwner!==session)clearIntake();intakeOwner=session;
    const schedule=frame.intake_schedule||{events:[]},signature=JSON.stringify(schedule);
    if(signature!==intakeSignature){intake.update(schedule);intakeSignature=signature;}
    syncIntake();onFrame(frame);
    const currentEnvironment={free:'studio',supine:'bed',upright:'floor'}[frame.mechanics?.body_environment?.kind];
    if(currentEnvironment)environment=currentEnvironment;
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
    status('Body running · live inputs act on the next tick.');
  }
  function pendingStatus(frame) {
    if(frame?.schema==='ihm.embodied-frame.v1'){initialized(frame);return;}
    if(frame?.error||frame?.closed||frame?.status==='error')throw Error(frame.error||'Body closed before initialization completed.');
    if(!['initializing','ready','closing'].includes(frame?.status))throw Error('Unexpected body startup state');
    initializing=true;syncIntake();status('Body initializing · waiting for the native resource slot.');
    panels.status('Body initializing · no live physiological frame yet.');
  }
  async function start() {
    if(disposed||resetting)return;
    if(faulted){status('Reset is required before this body can restart.');return;}
    if(creating)return creating;
    if(!session){
      creating=(async()=>{
        onPauseReplay();running=true;status('Initializing body…');
        const frame=await createBodyOwner(request,endpoint(),{environment:bodyEnvironment(environment,'embodied'),environment_selection:environmentSelection,...options});
        session=frame.id;
        if(!session)throw Error('Body startup returned no session identity');
        if(disposed){await request(endpoint()+'/'+session+'/close',{});return;}
        pendingStatus(frame);
      })();
      try{await creating;}catch(error){running=false;if(session)faulted=true;throw error;}finally{creating=null;}
    }else if(!initializing){running=true;panels.status('Body resumed · awaiting the next accepted native frame.');}
  }
  function pause(){running=false;panels.status('Paused · last accepted body state. Inputs take effect on resume.');}
  function reset() {
    if(resetTask)return resetTask;
    resetting=true;syncIntake();running=false;drag=null;gizmoDragging=false;arrow.visible=false;gizmo.detach();controls.enabled=true;
    resetTask=Promise.resolve().then(async()=>{
    try{
      if(creating){try{await creating;}catch(error){if(!session)throw error;}}
      while(pending||intakeTask)await new Promise(resolve=>setTimeout(resolve,10));
      running=false;
      if(session)await closeBodyOwner(request,endpoint()+'/'+session,async()=>{
        status('Closing body · waiting for native cleanup.');panels.status('Closing body · last accepted state remains shown until cleanup completes.');
        await new Promise(resolve=>setTimeout(resolve,500));
      });
      if(disposed)return;
      session=null;state=null;initializing=false;faulted=false;selectedId=null;panels.clear();spectrum.clear();intakeMass.clear();skinVoltage.clear();clearIntake();syncIntake();
      $('scene-clock').textContent='0.00 s';$('scene-force-value').textContent='0 N';
      onFrame(null);status('Body stopped.');
    }catch(e){status(e.message);throw e;}finally{resetting=false;resetTask=null;syncIntake();}
    });
    return resetTask;
  }
  $('scene-reset').onclick=()=>reset().catch(e=>status(e.message));
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
      catch(error){running=false;initializing=false;faulted=true;status(error.message);if(!disposed)panels.status('Body startup stopped: '+error.message);}
      finally{pending=false;}
      return;
    }
    if(!running||now-lastStep<20)return;
    pending=true;lastStep=now;
    try {
      const result=await advanceScene(request,endpoint()+'/'+session,bodyCommand('embodied',state.sequence,forceCommand(),panels.inputs,state));
      if(!disposed)accept(result.frame);
      if(result.recovered)throw Error(`${result.error.message}. Recorded state restored; resume when ready.`);
    }
    catch(e){running=false;drag=null;arrow.visible=false;controls.enabled=true;if(!disposed){panels.status('Stopped · last accepted body frame. '+e.message);if(e.message!==lastError)status(e.message);}lastError=e.message;}
    finally{pending=false;}
  }
  function unload(){if(session)navigator.sendBeacon(endpoint()+'/'+session+'/close',new Blob(['{}'],{type:'application/json'}));}
  window.addEventListener('pagehide',unload);
  return {update,reset,start,pause,scheduleIntakes,
    setEnvironment(name,selection={}){
      if(name===environment&&JSON.stringify(selection)===JSON.stringify(environmentSelection))return Promise.resolve();
      environment=name;environmentSelection=structuredClone(selection);
      return session||creating?reset():Promise.resolve();
    },
    setOptions(next){options={regional_skin:!!next.regional_skin,intake_mass:!!next.intake_mass};return session||creating?reset():Promise.resolve();},
    get running(){return running;},
    get started(){return !!session||!!creating;},
    get state(){return state;},
    get active(){return !!session||!!creating||pending;},
    dispose(){disposed=true;skinVoltage.dispose();intakeMass.dispose();spectrum.dispose();intake.dispose();running=false;unload();window.removeEventListener('pagehide',unload);gizmo.dispose();scene.remove(gizmo.getHelper());group.remove(environmentGroup,gizmoTarget,arrow);
      renderer.domElement.removeEventListener('pointerdown',down,true);renderer.domElement.removeEventListener('pointermove',move,true);renderer.domElement.removeEventListener('pointerup',up,true);renderer.domElement.removeEventListener('pointercancel',up,true);}};
}

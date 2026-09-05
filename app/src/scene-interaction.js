import * as THREE from 'three';
import { TransformControls } from 'three/addons/controls/TransformControls.js';
import './scene-interaction.css';
import {cursorSpring,rotateOffset,advanceScene} from './scene-forces.js';

export function mountSceneInteraction({scene,camera,renderer,controls,group,getObjects,getSelected,onSelect,onFrame,onPauseReplay}) {
  const mount=document.getElementById('scene-controls'),monitor=document.getElementById('scene-monitor');
  if(!mount||!monitor)throw Error('Scene interaction mounts missing');
  mount.innerHTML=`<label class="field-label" for="scene-environment">Environment</label>
    <select id="scene-environment"><option value="studio">Studio</option><option value="floor">Floor · supported stance</option><option value="bed">Bed · supported supine</option></select>
    <div class="scene-modes" role="group" aria-label="Cursor mode"><button data-scene-mode="select" aria-pressed="true">Select</button><button data-scene-mode="gimbal" aria-pressed="false">Gimbal</button><button data-scene-mode="force" aria-pressed="false">Force</button></div>
    <div class="scene-actions"><button id="scene-play">Start mechanics</button><button id="scene-reset">Reset</button></div>
    <p id="scene-status-note" class="muted" role="status">Select a part to inspect. Force mode applies a physical cursor spring.</p>`;
  monitor.innerHTML=`<div class="scene-values"><span>Time <b id="scene-clock">0.00 s</b></span><span>Applied force <b id="scene-force-value">0 N</b></span></div><p id="scene-target" class="muted">No force target</p><details class="evidence-fold"><summary>Mechanics and supports</summary><p id="scene-scope">Linked body translations and affine tissues; reference orientations constrained. Free objects have ground contact. Whole-body surface contact remains in development.</p></details>`;
  const $=id=>document.getElementById(id);
  const environmentGroup=new THREE.Group();environmentGroup.name='Interactive environment';group.add(environmentGroup);
  const objectMeshes=new Map(),ray=new THREE.Raycaster(),plane=new THREE.Plane(),cursor=new THREE.Vector3();
  const gizmoTarget=new THREE.Object3D();group.add(gizmoTarget);
  const gizmo=new TransformControls(camera,renderer.domElement);gizmo.setMode('translate');gizmo.setSize(.65);
  scene.add(gizmo.getHelper());
  const arrow=new THREE.ArrowHelper(new THREE.Vector3(0,1,0),new THREE.Vector3(),0,0xe9b979,.02,.009);group.add(arrow);arrow.visible=false;
  let mode='select',session=null,state=null,running=false,pending=false,creating=null,drag=null,gizmoDragging=false;
  let lastStep=0,disposed=false,environment='studio',selectedId=null,lastError='',resetting=false,resetTask=null;
  const status=text=>{if(!disposed)$('scene-status-note').textContent=text;};

  async function request(path,data) {
    const response=await fetch(path,data===undefined?{cache:'no-store'}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
    const payload=await response.json();if(!response.ok)throw Error(payload.error||'Scene request failed');return payload;
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
    if(environment!=='studio'){
      const bed=environment==='bed';
      const mesh=new THREE.Mesh(new THREE.BoxGeometry(bed?1.05:3,bed?2.1:.06,bed?.08:3),new THREE.MeshStandardMaterial({color:bed?0x718680:0x344646,roughness:.95}));
      if(bed)mesh.position.set(0,0,-.28);else mesh.position.set(0,-.99,0);
      environmentGroup.add(mesh);
    }
  }
  function accept(frame) {
    state=frame;onFrame(frame);
    const scope=frame.scope;
    if(scope)$('scene-scope').textContent=[scope.body_mechanics,scope.body_rotations,scope.body_gravity,scope.body_environment,
      scope.body_object_contact===false?'Body–object contact is not coupled.':null,
      scope.clothing_contact===false?'Clothing contact is not coupled to this scene.':null,
      scope.physiology_feedback===false?'This scene does not feed forces back into physiology.':null].filter(Boolean).join('. ');
    for(const item of frame.objects||[]){
      let mesh=objectMeshes.get(item.id);
      if(!mesh){mesh=new THREE.Mesh(new THREE.SphereGeometry(item.radius_m,24,16),new THREE.MeshStandardMaterial({color:0xd6a16e,roughness:.6}));mesh.userData.sceneObject=item.id;environmentGroup.add(mesh);objectMeshes.set(item.id,mesh);}
      mesh.position.fromArray(item.position_m);const r=item.rotation_matrix;mesh.setRotationFromMatrix(new THREE.Matrix4().set(r[0][0],r[0][1],r[0][2],0,r[1][0],r[1][1],r[1][2],0,r[2][0],r[2][1],r[2][2],0,0,0,0,1));
    }
    $('scene-clock').textContent=frame.time_s.toFixed(2)+' s';
    if(selectedId&&!gizmoDragging)gizmoTarget.position.copy(centroid(selectedId));
  }
  async function start() {
    if(disposed||resetting)return;
    if(creating)return creating;
    if(!session){
      creating=(async()=>{onPauseReplay();status('Initializing canonical mechanics…');const frame=await request('/api/scene/sessions',{environment});
        if(disposed){await request(`/api/scene/sessions/${frame.id}/close`,{});return;}
        session=frame.id;accept(frame);
        if(drag?.initialAnchor){
          if(drag.gimbal)drag.target.sub(drag.initialAnchor).add(centroid(drag.id));
          else drag.offset.copy(drag.initialAnchor).sub(centroid(drag.id));
          delete drag.initialAnchor;
        }
        status('Mechanics ready · supports and force records are retained.');})();
      try{await creating;}finally{creating=null;}
    }
    if(!disposed&&session){running=true;$('scene-play').textContent='Pause mechanics';}
  }
  function reset() {
    if(resetTask)return resetTask;
    resetting=true;running=false;drag=null;gizmoDragging=false;arrow.visible=false;gizmo.detach();controls.enabled=true;
    resetTask=Promise.resolve().then(async()=>{
    try{
      if(creating)await creating;
      while(pending)await new Promise(resolve=>setTimeout(resolve,10));
      running=false;
      if(session)await request(`/api/scene/sessions/${session}/close`,{});
      if(disposed)return;
      session=null;state=null;selectedId=null;environment=$('scene-environment').value;environmentMeshes();
      $('scene-play').textContent='Start mechanics';$('scene-clock').textContent='0.00 s';$('scene-force-value').textContent='0 N';
      onFrame(null);status('Scene reset. Start mechanics to apply forces.');
    }catch(e){status(e.message);throw e;}finally{resetting=false;resetTask=null;}
    });
    return resetTask;
  }
  $('scene-environment').onchange=()=>reset().catch(e=>status(e.message));$('scene-reset').onclick=()=>reset().catch(e=>status(e.message));
  $('scene-play').onclick=()=>{if(running){running=false;$('scene-play').textContent='Resume mechanics';}else start().catch(e=>status(e.message));};
  mount.querySelectorAll('[data-scene-mode]').forEach(button=>button.onclick=()=>{
    mode=button.dataset.sceneMode;mount.querySelectorAll('[data-scene-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b===button)));
    drag=null;arrow.visible=false;gizmo.detach();controls.enabled=true;
    if(mode==='gimbal'&&selectedId)gizmo.attach(gizmoTarget);
    status(mode==='select'?'Click a structure to inspect.':mode==='gimbal'?'Select a part; drag an axis to apply a cursor spring.':'Drag a body part or free object to apply force.');
  });
  gizmo.addEventListener('dragging-changed',event=>{
    gizmoDragging=event.value;controls.enabled=!event.value;
    if(event.value&&selectedId){drag={id:selectedId,offset:new THREE.Vector3(),target:gizmoTarget.position.clone(),initialAnchor:session?null:gizmoTarget.position.clone(),gimbal:true};start().catch(e=>{status(e.message);gizmoDragging=false;drag=null;});}
    else {drag=null;arrow.visible=false;}
  });
  gizmo.addEventListener('objectChange',()=>{if(gizmoDragging&&drag)drag.target.copy(gizmoTarget.position);});
  function down(event) {
    if(mode==='select'||event.button!==0||gizmoDragging||resetting)return;
    // Axis handles own their pointer events; only pick new anatomy away from an axis.
    if(mode==='gimbal'&&gizmo.axis)return;
    const hit=pick(event);if(!hit)return;selected(hit);
    if(mode==='gimbal')return;
    event.preventDefault();event.stopImmediatePropagation();controls.enabled=false;
    renderer.domElement.setPointerCapture(event.pointerId);
    plane.setFromNormalAndCoplanarPoint(camera.getWorldDirection(new THREE.Vector3()),hit.world);
    const rotation=state?.objects?.find(o=>o.id===hit.id)?.rotation_matrix;
    drag={id:hit.id,offset:new THREE.Vector3(...rotateOffset(hit.point.clone().sub(centroid(hit.id,hit.mesh)).toArray(),rotation,true)),target:hit.point.clone(),initialAnchor:session?null:hit.point.clone(),pointer:event.pointerId};
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
  function forceCommand() {
    if(!drag)return [];
    const rotation=state?.objects?.find(o=>o.id===drag.id)?.rotation_matrix;
    const point=centroid(drag.id).add(new THREE.Vector3(...rotateOffset(drag.offset.toArray(),rotation))),force=cursorSpring(point.toArray(),drag.target.toArray());
    const magnitude=Math.hypot(...force);$('scene-force-value').textContent=magnitude.toFixed(2)+' N';
    arrow.position.copy(point);arrow.visible=magnitude>1e-5;
    if(arrow.visible){arrow.setDirection(new THREE.Vector3(...force).normalize());arrow.setLength(Math.min(.3,magnitude*.015),.018,.009);}
    return [{id:drag.id,force_n:force,point_m:point.toArray()}];
  }
  async function update(now) {
    if(disposed||resetting||!running||pending||!session||now-lastStep<20)return;
    pending=true;lastStep=now;
    try {
      const result=await advanceScene(request,`/api/scene/sessions/${session}`,{seconds:.02,sequence:state.sequence,forces:forceCommand()});
      if(!disposed)accept(result.frame);
      if(result.recovered)throw Error(`${result.error.message}. Recorded state restored; resume when ready.`);
    }
    catch(e){running=false;drag=null;arrow.visible=false;controls.enabled=true;if(!disposed){$('scene-play').textContent='Resume mechanics';if(e.message!==lastError)status(e.message);}lastError=e.message;}
    finally{pending=false;}
  }
  function unload(){if(session)navigator.sendBeacon(`/api/scene/sessions/${session}/close`,new Blob(['{}'],{type:'application/json'}));}
  window.addEventListener('pagehide',unload);environmentMeshes();
  return {update,reset,createControls:()=>mount,get state(){return state;},get active(){return !!session;},
    dispose(){disposed=true;running=false;unload();window.removeEventListener('pagehide',unload);gizmo.dispose();scene.remove(gizmo.getHelper());group.remove(environmentGroup,gizmoTarget,arrow);
      renderer.domElement.removeEventListener('pointerdown',down,true);renderer.domElement.removeEventListener('pointermove',move,true);renderer.domElement.removeEventListener('pointerup',up,true);renderer.domElement.removeEventListener('pointercancel',up,true);}};
}

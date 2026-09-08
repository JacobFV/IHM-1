import * as THREE from "three";

// A pole-free orbit for a body that is examined from every side.
//
// OrbitControls is a turntable. It takes its orbit axis from `camera.up` once,
// in its constructor, and never looks at it again, so azimuth always runs about
// whatever `camera.up` was then -- here the body's superior axis. The polar
// angle is clamped to [0, PI] and then pushed off the poles by Spherical
// .makeSafe(), so the arc over the head, around the coronal plane, ends at a
// stop it cannot pass; and no drag rotates about the anterior-posterior axis at
// all. Walking around a standing body is free, looking at it from above and
// carrying on is not.
//
// Here a drag rotates the whole camera frame rigidly: dx about the camera's own
// up, dy about the camera's own right, both by an angle proportional to the
// pixels dragged and to nothing else. No axis is privileged, nothing is
// clamped, and there is no orientation in which a drag degenerates -- the arc
// over the head costs exactly what the arc around the waist costs.
//
// Rigid frame rotation is a trackball, and trackballs lean: two drags that do
// not commute leave the body's long axis tilted in the image. So while the
// pointer is idle the roll relaxes back toward `upReference`, the body's
// superior axis, projected into the image plane -- toward whichever of +/- that
// projection the current up is already nearer. That removes an accidental lean
// without undoing a deliberate trip over the top, which legitimately leaves the
// body inverted; a gimbal plane snap restores upright outright. The relaxation
// fades out as the view approaches the superior axis, where the projection is
// ill-conditioned and there is no lean to speak of.
//
// The API is the subset of OrbitControls this app uses -- `target`, `update()`,
// `enabled`, `enableDamping` -- plus the same pointer bindings, so
// scene-interaction's gizmo can still take the pointer by setting `enabled`.
const STATE = { NONE: 0, ROTATE: 1, DOLLY: 2, PAN: 3, TOUCH_ROTATE: 4, TOUCH_ZOOM_PAN: 5 };
const CHANGE = { type: "change" };

export function mountCameraOrbit(camera, domElement, options = {}) {
  const controls = new THREE.EventDispatcher();
  // The axis a settled view keeps vertical. It is the body's superior axis by
  // default, and whoever establishes an explicit view -- a plane snap, a scene
  // from the catalogue, which carries its own up -- sets it to that view's up so
  // relaxing never argues with the view it was given.
  controls.upReference = (options.upReference || new THREE.Vector3(0, 1, 0)).clone().normalize();

  controls.object = camera;
  controls.domElement = domElement;
  controls.target = new THREE.Vector3();
  controls.enabled = true;
  controls.enableRotate = true;
  controls.enableZoom = true;
  controls.enablePan = true;
  controls.enableDamping = false;
  controls.dampingFactor = 0.05;
  // One screen height of drag is one full turn, exactly as OrbitControls
  // measured it, so the fix does not quietly also change the gearing.
  controls.rotateSpeed = 1;
  controls.zoomSpeed = 1;
  controls.panSpeed = 1;
  controls.minDistance = 0;
  controls.maxDistance = Infinity;
  // Fraction of the remaining lean removed per frame while the pointer is idle,
  // and the sine of the angle to `upReference` below which no lean is removed.
  controls.rollRelaxFactor = 0.12;
  controls.rollRelaxFloor = 0.1;

  let state = STATE.NONE;
  const pending = new THREE.Vector2();   // radians still owed to the camera
  const panPending = new THREE.Vector3();
  let scalePending = 1;
  const pointers = [], positions = new Map();
  const dragStart = new THREE.Vector2();
  let touchStartDistance = 0;
  const touchPanStart = new THREE.Vector2();

  const offset = new THREE.Vector3();
  const back = new THREE.Vector3(), right = new THREE.Vector3(), frameUp = new THREE.Vector3();
  const spin = new THREE.Quaternion(), tilt = new THREE.Quaternion();
  const scratch = new THREE.Vector3(), wanted = new THREE.Vector3();
  const fallback = new THREE.Vector3();

  // The camera frame the renderer will actually use: three's lookAt builds it
  // from (position, target, up) exactly this way, so reading it back here and
  // rotating it keeps drags and pixels in step.
  function readFrame() {
    offset.copy(camera.position).sub(controls.target);
    const distance = offset.length();
    if (distance < 1e-9) return 0;
    back.copy(offset).divideScalar(distance);
    right.crossVectors(camera.up, back);
    if (right.lengthSq() < 1e-12) {
      // `up` collinear with the view axis leaves no right vector; any
      // perpendicular will do, and the next frame rebuilds from it.
      fallback.set(Math.abs(back.x) < 0.9 ? 1 : 0, Math.abs(back.x) < 0.9 ? 0 : 1, 0);
      right.crossVectors(fallback, back);
    }
    right.normalize();
    frameUp.crossVectors(back, right).normalize();
    return distance;
  }

  function height() {
    return domElement?.clientHeight || 1;
  }

  function rotate(dx, dy) {
    if (!controls.enableRotate) return;
    pending.x += (2 * Math.PI * dx * controls.rotateSpeed) / height();
    pending.y += (2 * Math.PI * dy * controls.rotateSpeed) / height();
  }
  // OrbitControls scales the orbit radius the frame it is asked to and never
  // damps the zoom; matching that keeps the wheel feeling identical.
  function dolly(factor) {
    if (!controls.enableZoom) return;
    scalePending *= factor;
  }
  const zoomScale = (delta) => Math.pow(0.95, controls.zoomSpeed * Math.abs(delta) * 0.01);
  function pan(dx, dy) {
    if (!controls.enablePan) return;
    const distance = readFrame();
    if (!distance) return;
    const reach = 2 * distance * Math.tan(((camera.fov || 50) / 2) * (Math.PI / 180)) * controls.panSpeed / height();
    panPending.addScaledVector(right, -dx * reach).addScaledVector(frameUp, dy * reach);
  }

  // ------------------------------------------------------------- pointers --
  function onPointerDown(event) {
    if (controls.enabled === false) return;
    if (pointers.length === 0) domElement.setPointerCapture?.(event.pointerId);
    pointers.push(event.pointerId);
    positions.set(event.pointerId, new THREE.Vector2(event.clientX, event.clientY));
    if (event.pointerType === "touch") {
      if (pointers.length === 1) { state = STATE.TOUCH_ROTATE; dragStart.set(event.clientX, event.clientY); }
      else if (pointers.length === 2) {
        state = STATE.TOUCH_ZOOM_PAN;
        const [a, b] = pointers.map((id) => positions.get(id));
        touchStartDistance = a.distanceTo(b);
        touchPanStart.set((a.x + b.x) / 2, (a.y + b.y) / 2);
      }
      return;
    }
    state = event.button === 1 ? STATE.DOLLY : event.button === 2 ? STATE.PAN : STATE.ROTATE;
    dragStart.set(event.clientX, event.clientY);
  }
  function onPointerMove(event) {
    if (controls.enabled === false || state === STATE.NONE) return;
    const previous = positions.get(event.pointerId);
    if (previous) previous.set(event.clientX, event.clientY);
    if (state === STATE.TOUCH_ZOOM_PAN) {
      if (pointers.length < 2) return;
      const [a, b] = pointers.map((id) => positions.get(id));
      const distance = a.distanceTo(b);
      if (touchStartDistance > 0 && distance > 0) dolly(touchStartDistance / distance);
      touchStartDistance = distance;
      const centre = new THREE.Vector2((a.x + b.x) / 2, (a.y + b.y) / 2);
      pan(centre.x - touchPanStart.x, centre.y - touchPanStart.y);
      touchPanStart.copy(centre);
      return;
    }
    const dx = event.clientX - dragStart.x, dy = event.clientY - dragStart.y;
    dragStart.set(event.clientX, event.clientY);
    if (state === STATE.ROTATE || state === STATE.TOUCH_ROTATE) rotate(dx, dy);
    else if (state === STATE.PAN) pan(dx, dy);
    else if (state === STATE.DOLLY && dy) dolly(dy > 0 ? 1 / zoomScale(dy) : zoomScale(dy));
  }
  function onPointerUp(event) {
    const index = pointers.indexOf(event.pointerId);
    if (index >= 0) pointers.splice(index, 1);
    positions.delete(event.pointerId);
    if (pointers.length === 0) { domElement.releasePointerCapture?.(event.pointerId); state = STATE.NONE; }
    else if (pointers.length === 1) { state = STATE.TOUCH_ROTATE; dragStart.copy(positions.get(pointers[0])); }
  }
  function onWheel(event) {
    if (controls.enabled === false || !controls.enableZoom) return;
    event.preventDefault();
    dolly(event.deltaY < 0 ? zoomScale(event.deltaY) : 1 / zoomScale(event.deltaY));
  }
  function onContextMenu(event) {
    if (controls.enabled !== false) event.preventDefault();
  }
  domElement.addEventListener("pointerdown", onPointerDown);
  domElement.addEventListener("pointermove", onPointerMove);
  domElement.addEventListener("pointerup", onPointerUp);
  domElement.addEventListener("pointercancel", onPointerUp);
  domElement.addEventListener("wheel", onWheel, { passive: false });
  domElement.addEventListener("contextmenu", onContextMenu);
  domElement.style.touchAction = "none";

  // ---------------------------------------------------------------- frame --
  controls.update = function update() {
    let distance = readFrame();
    if (!distance) return false;
    const share = controls.enableDamping ? controls.dampingFactor : 1;
    const yaw = pending.x * share, pitch = pending.y * share;
    pending.x -= yaw; pending.y -= pitch;
    if (Math.abs(pending.x) < 1e-9) pending.x = 0;
    if (Math.abs(pending.y) < 1e-9) pending.y = 0;
    let moved = false;

    if (yaw || pitch) {
      // Dragging right walks the camera left around the body, the same sense a
      // turntable gave: the body appears to follow the pointer.
      spin.setFromAxisAngle(frameUp, -yaw);
      tilt.setFromAxisAngle(right, -pitch);
      spin.multiply(tilt);
      offset.applyQuaternion(spin);
      camera.up.applyQuaternion(spin).normalize();
      camera.position.copy(controls.target).add(offset);
      moved = true;
      distance = readFrame();
    }
    if (scalePending !== 1) {
      const wantedDistance = THREE.MathUtils.clamp(distance * scalePending, controls.minDistance, controls.maxDistance);
      scalePending = 1;
      camera.position.copy(controls.target).addScaledVector(back, wantedDistance);
      moved = true;
      distance = readFrame();
    }
    if (panPending.lengthSq() > 0) {
      scratch.copy(panPending).multiplyScalar(controls.enableDamping ? controls.dampingFactor : 1);
      panPending.sub(scratch);
      if (panPending.lengthSq() < 1e-14) panPending.set(0, 0, 0);
      controls.target.add(scratch);
      camera.position.add(scratch);
      moved = true;
      distance = readFrame();
    }

    // Lean relaxation, only while the pointer is idle and nothing is owed.
    if (state === STATE.NONE && !pending.x && !pending.y) {
      const reference = controls.upReference;
      wanted.copy(reference).addScaledVector(back, -reference.dot(back));
      const lean = wanted.length();
      if (lean > controls.rollRelaxFloor) {
        wanted.divideScalar(lean);
        // Nearest upright, never a flip: a view reached by going over the head
        // is inverted on purpose and stays that way.
        if (wanted.dot(frameUp) < 0) wanted.negate();
        const angle = Math.atan2(scratch.crossVectors(frameUp, wanted).dot(back), frameUp.dot(wanted));
        if (Math.abs(angle) > 1e-5) {
          const fade = THREE.MathUtils.clamp((lean - controls.rollRelaxFloor) / 0.25, 0, 1);
          const step = Math.abs(angle) < 2e-3 ? angle : angle * controls.rollRelaxFactor * fade;
          if (step) {
            camera.up.applyQuaternion(spin.setFromAxisAngle(back, step)).normalize();
            readFrame();
            moved = true;
          }
        }
      }
    }

    camera.up.copy(frameUp);
    camera.lookAt(controls.target);
    if (moved) controls.dispatchEvent(CHANGE);
    return moved;
  };

  controls.dispose = function dispose() {
    domElement.removeEventListener("pointerdown", onPointerDown);
    domElement.removeEventListener("pointermove", onPointerMove);
    domElement.removeEventListener("pointerup", onPointerUp);
    domElement.removeEventListener("pointercancel", onPointerUp);
    domElement.removeEventListener("wheel", onWheel);
    domElement.removeEventListener("contextmenu", onContextMenu);
  };

  return controls;
}

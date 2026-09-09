import {canonicalEnvironmentObject} from '../../src/world-frame.js';
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';

test('real learned cortical controller advances native body and drags a bedroom ball',async({page},testInfo)=>{
 test.skip(process.env.RUN_NATIVE_BODY!=='1','Opt-in: initializes actual native physiology and IBM checkpoint');
 test.setTimeout(900000);page.setDefaultTimeout(30000);
 const frames=[],errors=[],commands=[],httpErrors=[];let owner=null;
 page.on('pageerror',error=>errors.push(error.message));
 page.on('request',request=>{if(request.method()==='POST'&&request.url().endsWith('/step'))commands.push(request.postDataJSON());});
 page.on('response',async response=>{
  if(!response.url().includes('/api/embodied/sessions'))return;
  try{const data=await response.json();if(!response.ok()||data.error){httpErrors.push(data.error||String(response.status()));console.log('Native browser HTTP error:',data.error);}if(data.id)owner=data.id;if(data.schema==='ihm.embodied-frame.v1'&&data.neural)frames.push(data);}catch{}
 });
 async function waitFrames(count,timeout){const end=Date.now()+timeout;while(frames.length<count){if(httpErrors.length)throw Error(httpErrors.join('; '));if(errors.length)throw Error(errors.join('; '));if(Date.now()>end)throw Error(`Timed out waiting for ${count} neural frames; received ${frames.length}`);await new Promise(resolve=>setTimeout(resolve,1000));}}
 await page.addInitScript(()=>localStorage.setItem('ihm.panes.v1',JSON.stringify({shown:['scene','live'],collapsed:[]})));
 await page.goto('/');console.log('Native browser: page loaded');
 await expect(page.locator('#scene-controller')).toHaveValue('regional',{timeout:60000});
 // The scene catalogue request shares the origin's connection pool with the
 // whole-body anatomy load (~2050 structures), so the environment tiles can
 // appear minutes after the controller select is ready. Wait for them
 // explicitly instead of letting the 30 s default decide.
 await expect(page.locator('#environment-tiles [data-tile="bedroom"]')).toBeVisible({timeout:300000});
 await page.locator('#environment-tiles [data-tile="bedroom"]').click();
 await page.locator('#insert-object').click();
 await page.locator('#environment-tiles [data-insert="ball-small"]').click();
 await page.locator('#scene-controller').selectOption('implicit_cortical_ankle');
 await page.locator('#scene-ablation').selectOption('no-cord');console.log('Native browser: learned cortical ankle no-cord bedroom+ball configured');
 try{
  await page.getByRole('button',{name:'Start body',exact:true}).click();
  await expect(page.locator('#scene-controller')).toBeDisabled();console.log('Native browser: owner creation requested');
  await waitFrames(1,600000);console.log('Native browser: first neural frame accepted');
  await waitFrames(3,120000);
  await page.getByRole('button',{name:'Pause body',exact:true}).click();
  const first=frames[0],last=frames.at(-1);
  expect(last.time_s).toBeGreaterThan(first.time_s);
  expect((last.controller||last.neural.controller).kind).toBe('implicit_cortical_ankle');
  expect(Object.keys(last.neural.arc_max).sort()).toEqual(['autogenic','reciprocal','renshaw','stretch']);
  expect(last.environment_state.objects.some(o=>o.kind==='cloth')).toBe(true);
  expect(last.environment_state.objects.some(o=>o.id.startsWith('ball-small'))).toBe(true);
  await expect(page.locator('.live-controller')).toContainText('Learned IBM cortical ankle control');
  await expect(page.locator('.live-controller')).toContainText('Stretch arc peak');
  const ball=canonicalEnvironmentObject(last.environment_state,last.environment_state.objects.find(o=>o.id.startsWith('ball-small')));
  await page.evaluate(position=>globalThis.__ihmCamera.setPose({target:position,position:position.map((v,i)=>v+(i===2?.7:0)),up:[0,1,0]}),ball.position_m);
  await page.locator('[data-scene-mode="force"]').click();
  const bounds=await page.locator('#scene').boundingBox(),x=bounds.x+bounds.width/2,y=bounds.y+bounds.height/2;
  await page.mouse.move(x,y);await page.mouse.down();await page.mouse.move(x+35,y-20,{steps:5});
  await expect.poll(()=>commands.some(c=>c.forces?.some(f=>f.id===ball.id&&Math.hypot(...f.force_n)>0)),{timeout:30000}).toBe(true);
  const previous=frames.length;
  await expect.poll(()=>frames.length,{timeout:60000}).toBeGreaterThan(previous);
  await page.mouse.up();
  await page.getByRole('button',{name:'Pause body',exact:true}).click();
  expect(errors).toEqual([]);
  await fs.writeFile(testInfo.outputPath('native-frames.json'),JSON.stringify({first,last},null,2));
  await page.screenshot({path:testInfo.outputPath('native-body.png')});
 }finally{
  await fs.writeFile(testInfo.outputPath('native-diagnostics.json'),JSON.stringify({errors,httpErrors,commands,frame_count:frames.length,first:frames[0],last:frames.at(-1)},null,2));
  if(!page.isClosed())await page.screenshot({path:testInfo.outputPath('native-final.png')}).catch(()=>{});
  if(owner)await page.request.post(`/api/embodied/sessions/${owner}/close`,{data:{}});
 }
});

test('real engineered balance advances upright native body without cortical claims',async({page},testInfo)=>{
 test.skip(process.env.RUN_NATIVE_STANCE!=='1','Opt-in: actual native engineered stance');
 test.setTimeout(900000);const frames=[],errors=[];let owner;
 page.on('pageerror',e=>errors.push(e.message));
 page.on('response',async response=>{if(!response.url().includes('/api/embodied/sessions'))return;try{const d=await response.json();if(d.error)errors.push(d.error);if(d.id)owner=d.id;if(d.schema==='ihm.embodied-frame.v1'&&d.neural)frames.push(d);}catch{}});
 await page.addInitScript(()=>localStorage.setItem('ihm.panes.v1',JSON.stringify({shown:['scene','live'],collapsed:[]})));
 await page.goto('/');await expect(page.locator('#scene-controller')).toHaveValue('regional',{timeout:60000});
 // Same catalogue/anatomy contention as the drag test above.
 await expect(page.locator('#environment-tiles [data-tile="floor"]')).toBeVisible({timeout:300000});
 await page.locator('#environment-tiles [data-tile="floor"]').click();
 await page.locator('#scene-controller').selectOption('engineering_stance');
 await expect(page.locator('#scene-ablation-label')).toBeHidden();
 await expect(page.locator('#scene-ankle-target-label')).toBeHidden();
 try{
  await page.getByRole('button',{name:'Start body',exact:true}).click();
  const until=Date.now()+600000;while(frames.length<3){if(errors.length)throw Error(errors.join('; '));if(Date.now()>until)throw Error(`Only ${frames.length} stance frames`);await new Promise(r=>setTimeout(r,1000));}
  await page.getByRole('button',{name:'Pause body',exact:true}).click();
  const last=frames.at(-1),controller={...last.neural.controller,...last.controller};
  expect(controller.kind).toBe('engineering_stance');expect(controller.motor_owner).toBe('native_local_discrete_lqr');
  expect(last.neural.lqr_stance).toBeTruthy();
  const surface=last.mechanics.surface_binding;
  expect(['ihm.segment-surface-binding.v1','ihm.continuous-surface-binding.v1']).toContain(surface.schema);
  expect(Object.keys(last.mechanics.surface_transforms)).toHaveLength(surface.segments.length);
  if(surface.weights_sha256){
   await expect.poll(()=>page.evaluate(async digest=>!!(await import('/src/surface-assets.js')).surfaceAssets.records.get(digest)?.asset,surface.weights_sha256),{timeout:60000}).toBe(true);
   await expect(page.locator('.live-controller')).toContainText('Shared graph-weighted native segment attachment');
  }else await expect(page.locator('.live-controller')).toContainText('Hard attachment to native segments');
  expect(last.time_s).toBeGreaterThan(frames[0].time_s);
  await expect(page.locator('.live-controller')).toContainText('Engineered balance');
  await expect(page.locator('.live-controller')).toContainText('Not used by this controller');expect(errors).toEqual([]);
 }finally{
  await fs.writeFile(testInfo.outputPath('stance-diagnostics.json'),JSON.stringify({errors,frame_count:frames.length,first:frames[0],last:frames.at(-1)},null,2));
  await page.screenshot({path:testInfo.outputPath('stance-final.png')});
  if(owner)await page.request.post(`/api/embodied/sessions/${owner}/close`,{data:{}});
 }
});

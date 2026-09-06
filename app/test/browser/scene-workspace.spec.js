import {test,expect} from '@playwright/test';

test('live scene owns spatial playback until reset and preserves the recorded run',async({page})=>{
  // A single triangle and two recorded samples keep this lifecycle test light.
  // Every API request is intercepted; no full-body assets or native jobs load.
  await page.route('**/api/**',async route=>{
    const path=new URL(route.request().url()).pathname;
    const fixtures={
      '/api/manifest':{models:[{id:'ihm-body',name:'Scene UI fixture',bounds:{min:[-.1,-.1,-.1],max:[.1,.1,.1]}}],structures:[{id:'fixture-body',name:'Fixture tissue',model_id:'ihm-body',system:'skeletal',kind:'mesh',default_visible:true}]},
      '/api/geometry/fixture-body':{positions:[-.08,0,0,.08,0,0,0,.08,0],indices:[0,1,2]},
      '/api/body':{name:'Explicit UI fixture',entity_count:1},
      '/api/body/trajectory':{centroids_m:{'fixture-body':[0,0,0]},frames:[{time_s:0,entities:{},physiology:{'ArterialPressure(mmHg)':90}},{time_s:1,entities:{},physiology:{'ArterialPressure(mmHg)':91}}]},
    };
    return route.fulfill(path in fixtures?{json:fixtures[path]}:{status:404,json:{error:'Unavailable in bounded UI fixture'}});
  });
  let sequence=0,closed=0,closing=false;
  const frame=()=>({id:'fixture-scene',schema:'ihm.embodied-frame.v1',model_id:'ihm-body',sequence,time_s:sequence*.02,
    entities:{'fixture-body':{centroid_m:[sequence*.001,0,0],translation_m:[sequence*.001,0,0],rotation_matrix:[[1,0,0],[0,1,0],[0,0,1]],deformation_gradient:[[1,0,0],[0,1,0],[0,0,1]]}},
    mechanics:{muscles:{bra_r:{activation:.1,excitation:.1,tendon_force_n:2,fiber_length_m:.1,sensor_basis:'UI fixture'}},body_environment:{kind:'supine',scope:'Controlled UI fixture'}} ,physiology:{values:{heart_rate_per_min:72+sequence},signal_metadata:{heart_rate_per_min:{label:'Heart rate',unit:'1/min',owner:'fixture'}}}});
  await page.route('**/api/embodied/sessions',r=>r.fulfill({json:r.request().method()==='POST'?{id:'fixture-scene',status:'initializing',closed:false}:{sessions:[]}}));
  await page.route('**/api/embodied/sessions/fixture-scene',r=>r.fulfill({json:closing?{id:'fixture-scene',closed:true,status:'closed'}:frame()}));
  await page.route('**/api/embodied/sessions/fixture-scene/step',r=>{sequence++;return r.fulfill({json:frame()});});
  await page.route('**/api/embodied/sessions/fixture-scene/close',r=>{closed++;closing=true;return r.fulfill({json:{closed:false,status:'closing'}});});
  const errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.goto('/');
  await expect(page.locator('#scene-play')).toBeVisible();
  await expect(page.locator('#play')).toBeEnabled({timeout:30000});
  const recorded=await page.locator('#trajectory-run').inputValue();
  await page.locator('#scene-play').click();
  await expect(page.locator('#scene-monitor')).toBeVisible();
  await expect.poll(()=>sequence).toBeGreaterThan(2);
  await expect(page.locator('#play')).toBeDisabled();
  await expect(page.locator('#time')).toBeDisabled();
  await expect(page.locator('#flow-legend')).toContainText('Unified live body');
  await expect(page.locator('#live-body-monitor')).toContainText('1 muscle effectors');
  await expect(page.locator('#live-signal-1-monitor')).toContainText('1/min');
  await expect(page.locator('#live-motor-monitor select')).toHaveValue('bra_r');
  await page.locator('#scene-reset').click();
  await expect.poll(()=>closed).toBe(1);
  await expect(page.locator('#play')).toBeEnabled();
  await expect(page.locator('#live-scene-signal-note')).toHaveAttribute('hidden','');
  await expect(page.locator('#trajectory-run')).toHaveValue(recorded);
  expect(errors).toEqual([]);
});

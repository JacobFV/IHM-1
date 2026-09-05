import {test,expect} from '@playwright/test';

function sparse(){return {model_id:'ihm-body',has_body_projection:false,configuration:{protocol:'meal',sample_interval_s:30},clock:{sample_interval_s:30},centroids_m:{},
  fields:{'Aorta.Glucose.concentration_mg_per_dl':{unit:'mg/dL'},'missing':{unit:'g'}},
  frames:[{time_s:0,entities:{},physiology:{},systemic_values:{'Aorta.Glucose.concentration_mg_per_dl':90,missing:null}},{time_s:30,entities:{},physiology:{},systemic_values:{'Aorta.Glucose.concentration_mg_per_dl':110,missing:null}}],
  spectra:{frequency_hz:[0,.01],sigma_per_s:[0,.1],variables:[{name:'Aorta.Glucose.concentration_mg_per_dl',removed_sample_mean:100,psd:{frequency_hz:[0,.01],psd:[0,2]},laplace_real:[[0,3],[0,6]],laplace_imag:[[0,4],[0,8]]}]},
  mechanism_edges:[{source:'stomach',target:'blood glucose',status:'native test record',implementation:'recorded transfer'}],limitations:['Test fixture, not human observations']};}

test('real native apnea record drives one dressed body, selected values and finite spectra',async({page,request})=>{
  const response=await request.get('/api/body/experiments/systemic-apnea');expect(response.ok()).toBeTruthy();
  const native=await response.json(),errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.goto('/');
  await expect(page.locator('#systemic-study option[value="apnea"]')).toHaveCount(1);
  await page.locator('#systemic-study').selectOption('apnea');
  await expect(page.locator('#systemic-monitor')).toBeVisible({timeout:120000});
  await expect(page.locator('#systemic-note')).toContainText('0.1 s recorded samples');
  await expect(page.locator('#systemic-note')).toContainText('computed respiratory body projection');
  await expect(page.locator('#trajectory-run')).toBeDisabled();
  await expect(page.locator('#trajectory-run')).toHaveValue('systemic:apnea');
  await expect(page.locator('#garment-shirt')).toBeChecked();await expect(page.locator('#garment-shorts')).toBeChecked();
  await page.locator('#variable').selectOption('arterial_co2_mmhg');
  await page.locator('#time').fill('900');await page.locator('#time').dispatchEvent('input');
  await expect(page.locator('#time-value')).toHaveText('90.000 s');
  await expect(page.locator('#systemic-selected')).toContainText(native.frames[900].systemic_values.arterial_co2_mmhg.toPrecision(5));
  await expect(page.locator('#chart .chart-values strong')).toHaveText(native.frames[900].systemic_values.arterial_co2_mmhg.toPrecision(5));
  await expect(page.locator('#flow-legend')).toContainText('thoracic skin');
  await page.locator('#tab-spectral').click();await expect(page.locator('#spectral-run')).toBeDisabled();
  await expect(page.locator('#spectral-run')).toHaveValue('systemic:apnea');
  await page.locator('#spectral-mode').selectOption('laplace');
  await expect(page.locator('#chart-note')).toContainText('Finite-horizon');
  await expect(page.locator('#chart .chart-values span')).toContainText('mmHg · s');
  await page.locator('#tab-phys').click();
  await expect(page.locator('#scene-status')).toHaveText('',{timeout:120000});
  await expect(page.locator('#time-value')).toHaveText('90.000 s');
  await page.screenshot({path:'test-results/systemic-apnea.png',fullPage:true});
  await page.locator('#systemic-study').selectOption('');
  await expect(page.locator('#systemic-monitor')).toBeHidden();
  await expect(page.locator('#trajectory-run')).toBeEnabled();
  await expect(page.locator('#chart-note')).toContainText('IHMGenericMale');
  expect(errors).toEqual([]);
});

test('late anatomy geometry preserves the selected systemic time',async({page})=>{
  let release;const gate=new Promise(resolve=>{release=resolve;});
  await page.route('**/api/manifest',r=>r.fulfill({json:{models:[{id:'ihm-body',name:'One body',bounds:{min:[-.5,-1,-.2],max:[.5,1,.2]}}],structures:[{id:'skin',name:'Skin',model_id:'ihm-body',system:'integumentary',kind:'mesh'}]}}));
  await page.route('**/api/geometry/skin',async r=>{await gate;await r.fulfill({json:{positions:[-.1,0,0,.1,0,0,0,.1,0],indices:[0,1,2]}});});
  await page.route('**/api/body/experiments/systemic',r=>r.fulfill({json:{runs:[{id:'meal',label:'Meal fixture',duration_s:30,sample_interval_s:30,has_body_projection:false}]}}));
  await page.route('**/api/body/experiments/systemic-meal',r=>r.fulfill({json:sparse()}));
  await page.goto('/');await page.locator('#systemic-study').selectOption('meal');
  await expect(page.locator('#systemic-monitor')).toBeVisible();
  await page.locator('#time').fill('1');await page.locator('#time').dispatchEvent('input');
  await expect(page.locator('#time-value')).toHaveText('30.000 s');release();
  await expect(page.locator('#scene-status')).toHaveText('');
  await expect(page.locator('#time-value')).toHaveText('30.000 s');
});

test('sparse meal data has explicit reference geometry and missing values remain unavailable',async({page})=>{
  const errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.route('**/api/body/experiments/systemic',r=>r.fulfill({json:{runs:[{id:'meal',label:'Meal fixture',duration_s:30,sample_interval_s:30,has_body_projection:false}]}}));
  await page.route('**/api/body/experiments/systemic-meal',r=>r.fulfill({json:sparse()}));
  await page.goto('/');await page.locator('#systemic-study').selectOption('meal');
  await expect(page.locator('#systemic-note')).toContainText('no body projection');
  await expect(page.locator('#systemic-speed')).toHaveValue('300');
  await page.locator('#time').fill('1');await page.locator('#time').dispatchEvent('input');
  await expect(page.locator('#flow-legend')).toContainText('reference anatomy');
  await expect(page.locator('#chart .chart-values strong')).toHaveText('110.00');
  await page.locator('#variable').selectOption('missing');
  await expect(page.locator('#systemic-selected')).toContainText('unavailable');
  await expect(page.locator('#chart .chart-values strong')).toHaveText('Unavailable');
  await page.locator('#tab-spectral').click();await page.locator('#spectral-mode').selectOption('laplace');await page.locator('#sigma').selectOption('1');
  await expect(page.locator('#chart .chart-values strong')).toHaveText('10.000');
  expect(errors).toEqual([]);
});

test('a late systemic response cannot replace a subsequently selected resting body',async({page})=>{
  let release;const gate=new Promise(resolve=>{release=resolve;});let reached=false;
  await page.route('**/api/body/experiments/systemic',r=>r.fulfill({json:{runs:[{id:'meal',label:'Late meal',duration_s:30,sample_interval_s:30,has_body_projection:false}]}}));
  await page.route('**/api/body/experiments/systemic-meal',async r=>{reached=true;await gate;await r.fulfill({json:sparse()});});
  await page.goto('/');await page.locator('#systemic-study').selectOption('meal');
  await expect.poll(()=>reached).toBe(true);
  await expect(page.locator('#play')).toBeDisabled();
  await page.locator('#systemic-study').selectOption('');release();
  await expect(page.locator('#systemic-monitor')).toBeHidden();
  await expect(page.locator('#chart-note')).toContainText('IHMGenericMale');
  await expect(page.locator('#systemic-study')).toHaveValue('');
});

test('late resting trajectory cannot reset an active regional electrical clock',async({page})=>{
  let release;const gate=new Promise(resolve=>{release=resolve;});let completed=false;
  await page.route('**/api/body/trajectory*',async r=>{await gate;await r.fulfill({json:{centroids_m:{},frames:[{time_s:0,entities:{},physiology:{}},{time_s:1,entities:{},physiology:{}}]}});completed=true;});
  await page.goto('/');await expect(page.locator('#model')).toHaveValue('ihm-body');
  await page.locator('#regional-study').selectOption('skin-electric');
  await expect(page.locator('#details h2')).toHaveText('Non-neural skin electricity');
  await page.locator('#time').fill('40');await page.locator('#time').dispatchEvent('input');
  await expect(page.locator('#flow-legend')).toContainText('92.74 V/m');release();
  await expect.poll(()=>completed).toBe(true);
  await expect(page.locator('#time')).toHaveValue('40');
  await expect(page.locator('#flow-legend')).toContainText('92.74 V/m');
});

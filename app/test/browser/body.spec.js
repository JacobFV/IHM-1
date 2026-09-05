import {test,expect} from '@playwright/test';
const trajectory = {centroids_m:{'body-bp3d-FJ3259':[.1,-.3,0]},frames:[
  {time_s:0,entities:{},physiology:{}},
  {time_s:1,entities:{'body-bp3d-FJ3259':{translation_m:[.001,0,0]}},physiology:{'MeanArterialPressure(mmHg)':90}}
]};
test('one canonical body opens by default with explicit evidence and computed playback',async({page})=>{
  const errors=[];page.on('pageerror',error=>errors.push(error.message));
  await page.route('**/api/body',route=>route.fulfill({json:{name:'IHM',entity_count:2408,status:'assembled'}}));
  await page.route('**/api/body/trajectory',route=>route.fulfill({json:trajectory}));
  await page.goto('/');
  await expect(page.locator('#model')).toHaveValue('ihm-body');
  await expect(page.locator('#source-inspection')).not.toHaveAttribute('open');
  await expect(page.locator('#patient')).toHaveValue('IHMGenericMale');
  await expect(page.locator('#patient')).toBeDisabled();
  await expect(page.locator('#scene-status')).toHaveText('',{timeout:120000});
  await page.locator('#search').fill('superior lobe of left lung');
  await page.locator('#structures button').first().click();
  await expect(page.locator('#details')).toContainText('Biological uncertainty');
  await expect(page.locator('#details')).toContainText('Registration & synthesis');
  await expect(page.locator('#details')).toContainText('Display & numerical precision');
  await expect(page.locator('#details')).toContainText('4.67 mm RMS');
  await expect(page.locator('#play')).toBeEnabled();
  await page.locator('#time').fill('1');await page.locator('#time').dispatchEvent('input');
  await expect(page.locator('#time-value')).toHaveText('1.000 s');
  await expect(page.locator('#flow-legend')).toContainText('MAP 90.0 mmHg');
  await page.locator('#source-inspection summary').click();
  await page.locator('#model').selectOption('bodyparts3d');
  await expect(page.locator('#patient')).toBeEnabled();
  await page.locator('#canonical-body').click();
  await expect(page.locator('#model')).toHaveValue('ihm-body');
  await expect(page.locator('#source-inspection')).not.toHaveAttribute('open');
  expect(errors).toEqual([]);
});
test('unavailable body trajectory never substitutes animated anatomy',async({page})=>{
  await page.route('**/api/body/trajectory',route=>route.fulfill({status:404,json:{error:'Body trajectory not computed'}}));
  await page.goto('/');
  await expect(page.locator('#model')).toHaveValue('ihm-body');
  await expect(page.locator('#time-value')).toHaveText('Body trajectory unavailable');
  await expect(page.locator('#play')).toBeDisabled();
  await expect(page.locator('#time')).toBeDisabled();
});
test('real canonical trajectory drives the body clock and physiology without extrapolation',async({page,request})=>{
  const response=await request.get('/api/body/trajectory');
  expect(response.ok()).toBeTruthy();
  const data=await response.json();
  expect(data.frames.length).toBeGreaterThan(1);
  const errors=[];page.on('pageerror',error=>errors.push(error.message));
  await page.goto('/');
  await expect(page.locator('#model')).toHaveValue('ihm-body');
  await expect(page.locator('#scene-status')).toHaveText('',{timeout:120000});
  await expect(page.locator('#play')).toBeEnabled();
  await expect(page.locator('#trajectory-run')).toHaveValue('body');
  await expect(page.locator('#chart-note')).toContainText('IHMGenericMale');
  const index=Math.min(50,data.frames.length-1);
  await page.locator('#time').fill(String(index));await page.locator('#time').dispatchEvent('input');
  await expect(page.locator('#time-value')).toHaveText(`${data.frames[index].time_s.toFixed(3)} s`);
  await expect(page.locator('#flow-legend')).toContainText(`${Object.keys(data.frames[index].entities).length} tissue transforms`);
  await page.locator('#search').fill('superior lobe of left lung');
  await page.locator('#structures button').first().click();
  await expect(page.locator('#details')).toContainText('registered geometry');
  await page.locator('#details summary').click();
  await expect(page.locator('#details')).toContainText('Shared named bone');
  await page.locator('#search').fill('');
  await expect(page.locator('#scene-status')).toHaveText('',{timeout:120000});
  await page.screenshot({path:'test-results/canonical-body.png',fullPage:true});
  expect(errors).toEqual([]);
});

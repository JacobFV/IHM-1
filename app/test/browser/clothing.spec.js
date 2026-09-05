import {test,expect} from '@playwright/test';

test('canonical human wears independent opaque clothing by default and supports anatomy inspection',async({page})=>{
  const errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.goto('/');
  await expect(page.locator('#garment-shirt')).toBeChecked();
  await expect(page.locator('#garment-shorts')).toBeChecked();
  await expect(page.locator('#clothing-status')).toContainText('Fitted to canonical skin',{timeout:120000});
  await expect(page.locator('#scene-status')).toHaveText('',{timeout:120000});
  await page.screenshot({path:'test-results/clothed-human.png',fullPage:true});
  const dressed=await page.locator('#viewport canvas').screenshot();
  await page.locator('#garment-shirt').uncheck();
  await page.locator('#garment-shorts').uncheck();
  await expect(page.locator('#clothing-status')).toContainText('0 garments shown');
  const anatomy=await page.locator('#viewport canvas').screenshot();
  expect(anatomy.equals(dressed)).toBe(false);
  await page.getByLabel('Anatomy view',{exact:true}).selectOption('core');
  await expect(page.locator('#garment-shirt')).not.toBeChecked();
  await page.locator('#garment-shorts').check();
  await expect(page.locator('#clothing-status')).toContainText('1 garment shown');
  await page.locator('#source-inspection summary').click();
  await page.locator('#model').selectOption('bodyparts3d');
  await expect(page.locator('#clothing-components')).toBeHidden();
  await page.locator('#canonical-body').click();
  await expect(page.locator('#clothing-components')).toBeVisible();
  await expect(page.locator('#garment-shirt')).not.toBeChecked();
  await expect(page.locator('#garment-shorts')).toBeChecked();
  expect(errors).toEqual([]);
});

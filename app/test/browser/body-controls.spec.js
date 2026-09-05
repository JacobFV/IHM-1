import { test, expect } from '@playwright/test';

test('body interventions and resizable independent panels', async ({ page }) => {
  await page.goto('/');
  await expect(page.locator('#body-protocol')).toBeVisible();
  await expect(page.locator('#chest-compliance')).toBeVisible();
  const handle=page.getByRole('separator',{name:'Resize anatomy panel'});
  await expect(handle).toBeVisible();
  const before=await page.locator('#library-panel').boundingBox();
  await handle.focus();
  await page.keyboard.press('ArrowRight');
  const after=await page.locator('#library-panel').boundingBox();
  expect(after.width).toBeGreaterThan(before.width);
  expect(await page.evaluate(()=>document.documentElement.scrollHeight<=innerHeight)).toBe(true);
  await page.getByRole('button',{name:'Anatomy panel',exact:true}).click();
  await expect(handle).toBeHidden();
});

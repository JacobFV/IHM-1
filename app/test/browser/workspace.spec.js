import { test, expect } from '@playwright/test';

async function expectNoPageOverflow(page) {
  expect(await page.evaluate(() => ({
    width: document.documentElement.scrollWidth <= innerWidth,
    height: document.documentElement.scrollHeight <= innerHeight,
    top: scrollY,
  }))).toEqual({ width: true, height: true, top: 0 });
}

test('the scene expands when panels close and focus restores panel choices and inputs', async ({ page }) => {
  await page.goto('/');
  await expect(page.locator('#viewport canvas')).toBeVisible();
  await expectNoPageOverflow(page);
  const before = await page.locator('#viewport').boundingBox();
  await page.locator('#search').fill('lung');
  await page.locator('#duration').fill('17');
  await page.getByRole('button', { name: 'Anatomy panel', exact: true }).click();
  await expect(page.locator('#library-panel')).toBeHidden();
  await expect(page.getByRole('button', { name: 'Anatomy panel', exact: true })).toHaveAttribute('aria-expanded', 'false');
  expect((await page.locator('#viewport').boundingBox()).width).toBeGreaterThan(before.width);
  await page.getByRole('button', { name: 'Focus mode', exact: true }).click();
  await expect(page.locator('#inspector-panel')).toBeHidden();
  await expect(page.locator('#signals-panel')).toBeHidden();
  const focused = await page.locator('#viewport').boundingBox();
  expect(focused.width).toBeGreaterThan(before.width);
  expect(focused.height).toBeGreaterThan(before.height);
  await page.getByRole('button', { name: 'Focus mode', exact: true }).click();
  await expect(page.locator('#library-panel')).toBeHidden();
  await expect(page.locator('#inspector-panel')).toBeVisible();
  await expect(page.locator('#signals-panel')).toBeVisible();
  await expect(page.locator('#duration')).toHaveValue('17');
  await page.getByRole('button', { name: 'Anatomy panel', exact: true }).click();
  await expect(page.locator('#search')).toHaveValue('lung');
  await page.getByRole('button', { name: 'Signals panel', exact: true }).click();
  await page.reload();
  await expect(page.locator('#signals-panel')).toBeHidden();
  await expect(page.locator('#library-panel')).toBeVisible();
  await expectNoPageOverflow(page);
});

test('each panel scrolls within the fixed workspace and hiding it retains its scroll position', async ({ page }) => {
  await page.setViewportSize({ width: 1280, height: 650 });
  await page.goto('/');
  await expect(page.locator('#systems input').first()).toBeVisible();
  await page.locator('#inspector-panel .evidence-fold').evaluateAll(nodes => nodes.forEach(node => node.open = true));
  for (const id of ['library-panel', 'inspector-panel', 'signals-panel']) {
    const panel = page.locator(`#${id}`);
    const original = await panel.evaluate(el => { el.scrollTop = el.scrollHeight; return el.scrollTop; });
    expect(original).toBeGreaterThan(0);
    await expectNoPageOverflow(page);
    if (id === 'library-panel') {
      await page.getByRole('button', { name: 'Anatomy panel', exact: true }).click();
      await page.getByRole('button', { name: 'Anatomy panel', exact: true }).click();
      expect(await panel.evaluate(el => el.scrollTop)).toBe(original);
    }
  }
  const box = await page.locator('#viewport').boundingBox();
  expect(box.height).toBeGreaterThan(260);
});

test('narrow screens keep a useful scene and keyboard-operable drawers without page scrolling', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('/');
  await expect(page.locator('#viewport canvas')).toBeVisible();
  await expect(page.locator('#library-panel')).toBeHidden();
  await expect(page.locator('#inspector-panel')).toBeHidden();
  const anatomyToggle = page.getByRole('button', { name: 'Anatomy panel', exact: true });
  await anatomyToggle.focus();
  await page.keyboard.press('Enter');
  await expect(page.locator('#library-panel')).toBeVisible();
  await page.locator('#search').fill('lung');
  await page.keyboard.press('Escape');
  await expect(page.locator('#library-panel')).toBeHidden();
  await expect(anatomyToggle).toBeFocused();
  await page.getByRole('button', { name: 'Inspector panel', exact: true }).click();
  await expect(page.locator('#inspector-panel')).toBeVisible();
  await anatomyToggle.click();
  await expect(page.locator('#library-panel')).toBeVisible();
  await expect(page.locator('#inspector-panel')).toBeHidden();
  await expect(page.locator('#search')).toHaveValue('lung');
  await page.getByRole('button', { name: 'Focus mode', exact: true }).click();
  await expect(page.locator('#library-panel')).toBeHidden();
  expect((await page.locator('#viewport').boundingBox()).height).toBeGreaterThan(500);
  await expectNoPageOverflow(page);
  await page.setViewportSize({ width: 320, height: 568 });
  await expectNoPageOverflow(page);
  await expect(page.getByRole('button', { name: 'Focus mode', exact: true })).toBeVisible();
});

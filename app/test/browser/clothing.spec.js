import { test, expect } from "@playwright/test";

test("clothing tiles come from the garment catalog and change the render", async ({ page, request }) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  const catalog = await (await request.get("/api/clothing")).json();
  expect(catalog.garments.length).toBeGreaterThan(0);
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  await expect(page.locator("#clothing-tiles .tile span")).toHaveText(catalog.garments.map((g) => g.label));
  for (const garment of catalog.garments)
    await expect(page.locator(`#clothing-tiles [data-tile="${garment.id}"]`)).toHaveAttribute("aria-pressed", "true");
  const dressed = await page.locator("#scene").screenshot();
  for (const garment of catalog.garments)
    await page.locator(`#clothing-tiles [data-tile="${garment.id}"]`).click();
  for (const garment of catalog.garments)
    await expect(page.locator(`#clothing-tiles [data-tile="${garment.id}"]`)).toHaveAttribute("aria-pressed", "false");
  await page.waitForTimeout(500);
  expect((await page.locator("#scene").screenshot()).equals(dressed)).toBe(false);
  expect(errors).toEqual([]);
});

test("garments sharing a catalog slot are mutually exclusive", async ({ page }) => {
  await page.route("**/api/clothing", (route) => route.fulfill({ json: {
    schema: "ihm.clothing-catalog.v1",
    garments: [
      { id: "shirt", label: "Sleeveless shirt", slot: "torso-base", thumbnail_url: null },
      { id: "shorts", label: "Shorts", slot: "legs", thumbnail_url: null },
    ] } }));
  await page.goto("/");
  await expect(page.locator("#clothing-tiles .tile")).toHaveCount(2);
  await page.locator('#clothing-tiles [data-tile="shirt"]').click();
  await expect(page.locator('#clothing-tiles [data-tile="shirt"]')).toHaveAttribute("aria-pressed", "false");
  // Different slots stay independent.
  await expect(page.locator('#clothing-tiles [data-tile="shorts"]')).toHaveAttribute("aria-pressed", "true");
});

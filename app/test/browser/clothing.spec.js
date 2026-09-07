import { test, expect } from "@playwright/test";

test("clothing tiles come from the wardrobe and change the render", async ({ page, request }) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  const catalog = await (await request.get("/api/clothing")).json();
  expect(catalog.garments.length).toBeGreaterThan(2);
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  await expect(page.locator("#clothing-tiles .tile span")).toHaveText(catalog.garments.map((g) => g.label));

  // Slots contend, so an opening outfit wears some garments and not all of them.
  const worn = page.locator('#clothing-tiles [aria-pressed="true"]');
  const count = await worn.count();
  expect(count).toBeGreaterThan(0);
  expect(count).toBeLessThan(catalog.garments.length);

  const dressed = await page.locator("#scene").screenshot();
  for (const id of await worn.evaluateAll((els) => els.map((e) => e.dataset.tile)))
    await page.locator(`#clothing-tiles [data-tile="${id}"]`).click();
  await expect(page.locator('#clothing-tiles [aria-pressed="true"]')).toHaveCount(0);
  await page.waitForTimeout(500);
  expect((await page.locator("#scene").screenshot()).equals(dressed)).toBe(false);
  expect(errors).toEqual([]);
});

test("a garment excludes anything sharing any of its slots", async ({ page, request }) => {
  const catalog = await (await request.get("/api/clothing")).json();
  // Use the real wardrobe: a dress spans torso and legs, so it must clear both.
  const spanning = catalog.garments.find((g) => (g.slots || []).length > 1);
  test.skip(!spanning, "no multi-slot garment in the wardrobe");
  const [a, b] = spanning.slots;
  const inA = catalog.garments.find((g) => g.id !== spanning.id && (g.slots || []).includes(a) && (g.slots || []).length === 1);
  const inB = catalog.garments.find((g) => g.id !== spanning.id && (g.slots || []).includes(b) && (g.slots || []).length === 1);
  const other = catalog.garments.find((g) => !(g.slots || []).some((s) => spanning.slots.includes(s)));
  test.skip(!inA || !inB || !other, "wardrobe lacks the single-slot garments this needs");

  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  const tile = (id) => page.locator(`#clothing-tiles [data-tile="${id}"]`);
  const pressed = async (id) => await tile(id).getAttribute("aria-pressed");

  for (const g of [inA, inB, other]) if ((await pressed(g.id)) !== "true") await tile(g.id).click();
  expect(await pressed(inA.id)).toBe("true");
  expect(await pressed(inB.id)).toBe("true");

  await tile(spanning.id).click();
  expect(await pressed(spanning.id)).toBe("true");
  expect(await pressed(inA.id)).toBe("false");
  expect(await pressed(inB.id)).toBe("false");
  expect(await pressed(other.id)).toBe("true");
});

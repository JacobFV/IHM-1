import { test, expect } from "@playwright/test";

const trajectory = {
  centroids_m: { "body-bp3d-FJ3259": [0.1, -0.3, 0] },
  frames: [
    { time_s: 0, entities: {}, physiology: {} },
    { time_s: 1, entities: { "body-bp3d-FJ3259": { translation_m: [0.001, 0, 0] } }, physiology: {} },
  ],
};

// The distinction the redesign exists to make: the transport is the
// simulation's until a recording is explicitly opened, and while one is open
// the header names it and the transport scrubs it. No shared button, no mode
// toggle, and nothing carried by a tooltip.
test("the transport belongs to the simulation until a recording is opened", async ({ page }) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.route("**/api/body/trajectory*", (route) => route.fulfill({ json: trajectory }));
  await page.goto("/");
  await expect(page.locator("#transport")).toHaveAttribute("data-owner", "simulation");
  await expect(page.locator("#transport-label")).toHaveText("Simulation");
  await expect(page.locator("#play")).toHaveAttribute("title", /Run the simulation/);
  // The old control's title said "This does not run the simulation". Nothing
  // in this app says that any more, because nothing needs to.
  await expect(page.locator("#play")).not.toHaveAttribute("title", /does not run/);
  await expect(page.locator("#recording-header")).toBeHidden();
  await expect(page.locator("#time")).toBeHidden();
  await expect(page.locator("#speed")).toBeHidden();
  await expect(page.locator("#time-value")).toHaveText(/No recording is open/);

  await page.locator("#open-recording").click({ force: true });
  await page.locator("#recording-menu button").first().click({ force: true });
  await expect(page.locator("#transport")).toHaveAttribute("data-owner", "recording", { timeout: 240000 });
  await expect(page.locator("#transport-label")).toHaveText("Recording");
  await expect(page.locator("#recording-header")).toBeVisible();
  await expect(page.locator("#recording-name")).toHaveText("Canonical trajectory");
  await expect(page.locator("#time")).toBeVisible();
  await expect(page.locator("#speed")).toBeVisible();
  await page.locator("#time").fill("1");
  await page.locator("#time").dispatchEvent("input");
  await expect(page.locator("#time-value")).toHaveText("Canonical trajectory · recorded · 1.000 s");
  await expect(page.locator("#recording-detail")).toHaveText(/frame 2 \/ 2/);

  // Closing hands the transport straight back to the live body.
  await page.locator("#recording-close").click({ force: true });
  await expect(page.locator("#transport")).toHaveAttribute("data-owner", "simulation");
  await expect(page.locator("#recording-header")).toBeHidden();
  await expect(page.locator("#time")).toBeHidden();
  expect(errors).toEqual([]);
});

test("an unavailable recording says why rather than looking broken", async ({ page }) => {
  await page.route("**/api/body/trajectory*", (route) =>
    route.fulfill({ status: 409, json: { error: "Canonical trajectory is stale; rematerialize the native run" } }));
  await page.goto("/");
  // A stale trajectory is not a broken app: the simulation's own transport is
  // untouched by it, and the recording that will not load is never offered.
  await expect(page.locator("#transport")).toHaveAttribute("data-owner", "simulation", { timeout: 240000 });
  await expect(page.locator("#play")).toBeEnabled();
  await page.locator("#open-recording").click({ force: true });
  await expect(page.locator("#recording-menu")).not.toContainText("Reading what has been computed");
  await expect(page.locator("#recording-menu")).not.toContainText("Canonical trajectory");
  const entries = page.locator("#recording-menu button");
  if (await entries.count()) {
    // Opening one that will not load still names it, and says why instead of
    // reciting a frame count the file does not deliver.
    await entries.first().click({ force: true });
    await expect(page.locator("#recording-header")).toHaveAttribute("data-failed", "true");
    await expect(page.locator("#recording-detail")).toHaveText(/stale/);
    await expect(page.locator("#play")).toBeDisabled();
    await expect(page.locator("#play")).toHaveAttribute("title", /stale/);
  } else {
    await expect(page.locator("#recording-menu")).toContainText(/No recording is available|No computed recording/);
  }
});

test("the ring draws a dead path dead and the load-bearing one thick", async ({ page }) => {
  await page.goto("/");
  const items = page.locator(".ring-item");
  await expect(items.first()).toBeVisible({ timeout: 240000 });
  // With no body running the basis is stored measurement, and it says so.
  await expect(page.locator("#ring-basis")).toHaveAttribute("data-live", "false");
  await expect(page.locator("#ring-basis")).toContainText(/Measured, not live/);

  await page.locator('#ring-systems button[data-system="cord"]').click({ force: true });
  const widths = await page.locator("#ring-panel .ring-rank .ring-gauge line")
    .evaluateAll((nodes) => nodes.map((n) => Number(n.getAttribute("stroke-width"))));
  expect(widths.length).toBeGreaterThan(4);
  // The severed-cortex result is the one thick arrow; everything below it is a
  // hairline, because everything below it was measured to be one.
  expect(Math.max(...widths)).toBeGreaterThan(10);
  expect(widths.filter((w) => w <= 1).length).toBeGreaterThan(3);
  await expect(page.locator("#ring-panel")).toContainText("100.0% · measured");
  await expect(page.locator("#ring-panel")).toContainText("not measured");
});

test("annotations stand on their own projected points and groups toggle", async ({ page }) => {
  await page.goto("/");
  const items = page.locator(".ring-item:not([hidden])");
  await expect(items.first()).toBeVisible({ timeout: 240000 });

  // Every label carries the point it is about: either its own edge's, or its
  // system's anchor with the fallback declared rather than implied.
  const anchors = await page.locator(".ring-item").evaluateAll(
    (nodes) => nodes.map((n) => n.dataset.anchor));
  expect(anchors.length).toBeGreaterThan(0);
  expect(anchors.every((a) => a === "own" || a === "system")).toBe(true);

  // Groups are a multi-select now: turning them all on shows more than one.
  const single = await items.count();
  await page.locator("#ring-all").click({ force: true });
  await expect
    .poll(async () => items.count(), { timeout: 20000 })
    .toBeGreaterThan(single);

  // A label is placed where its point projects, so orbiting the body moves it.
  // A label that leaves the view is hidden and keeps its last transform, so the
  // question is whether the SET of placements moved, not any one of them.
  const places = () => page.locator(".ring-item:not([hidden])")
    .evaluateAll((nodes) => nodes.map((n) => n.style.transform).join("|"));
  const before = await places();
  const view = page.locator("canvas").first();
  const bounds = await view.boundingBox();
  await page.mouse.move(bounds.x + bounds.width / 2, bounds.y + bounds.height / 3);
  await page.mouse.down();
  await page.mouse.move(bounds.x + bounds.width / 2 + 260, bounds.y + bounds.height / 3, { steps: 16 });
  await page.mouse.up();
  await expect.poll(places, { timeout: 15000 }).not.toBe(before);

  await page.screenshot({ path: "test-results/annotations.png", fullPage: false });

  // Clicking a shown-and-focused group hides it.
  const cord = page.locator('#ring-systems button[data-system="cord"]');
  await cord.click({ force: true });                      // focus it
  await expect(cord).toHaveAttribute("aria-pressed", "true");
  await cord.click({ force: true });                      // hide it
  await expect(cord).toHaveAttribute("aria-pressed", "false");
  await expect(page.locator('.ring-item[data-group="cord"]:not([hidden])')).toHaveCount(0);
});

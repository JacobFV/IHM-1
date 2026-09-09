import { test, expect } from "@playwright/test";

const trajectory = {
  centroids_m: { "body-bp3d-FJ3259": [0.1, -0.3, 0] },
  frames: [
    { time_s: 0, entities: {}, physiology: {} },
    { time_s: 1, entities: { "body-bp3d-FJ3259": { translation_m: [0.001, 0, 0] } }, physiology: {} },
  ],
};

test("a computed trajectory drives the transport and the playback pane", async ({ page }) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.route("**/api/body/trajectory*", (route) => route.fulfill({ json: trajectory }));
  await page.goto("/");
  await expect(page.locator("#play")).toBeEnabled({ timeout: 240000 });
  await page.locator("#time").fill("1");
  await page.locator("#time").dispatchEvent("input");
  await expect(page.locator("#time-value")).toHaveText("Recorded \u00b7 1.000 s");
  // The bottom-centre control names its tenant, so it cannot be read as the
  // simulation's play button.
  await expect(page.locator("#transport-label")).toHaveText("Recording");
  await expect(page.locator("#transport")).toHaveAttribute("data-empty", "false");
  await expect(page.locator("#live-readout")).toBeHidden();
  expect(errors).toEqual([]);
});

test("an unavailable trajectory never substitutes animated anatomy", async ({ page }) => {
  await page.route("**/api/body/trajectory*", (route) =>
    route.fulfill({ status: 409, json: { error: "Canonical trajectory is stale; rematerialize the native run" } }));
  await page.goto("/");
  await expect(page.locator("#time-value")).toHaveText("Canonical trajectory is stale; rematerialize the native run");
  await expect(page.locator("#play")).toBeDisabled();
  await expect(page.locator("#time")).toBeDisabled();
  // A dead play button must say why it is dead rather than look broken.
  await expect(page.locator("#transport-label")).toHaveText("No recording");
  await expect(page.locator("#transport")).toHaveAttribute("data-empty", "true");
  await expect(page.locator("#play")).toHaveAttribute("title", /No recording loaded/);
});

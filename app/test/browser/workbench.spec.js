import { test, expect } from "@playwright/test";
test("real source anatomy renders and research controls remain connected", async ({
  page,
}) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto("/");
  await expect(page.locator("#model option")).not.toHaveText([
    "Loading source families…",
  ]);
  await expect(page.locator("#viewport canvas")).toBeVisible();
  await expect(page.locator("#render-count")).not.toHaveText(
    "0 structures visible",
    { timeout: 120000 },
  );
  await expect(page.locator("#scene-status")).toHaveText("", {
    timeout: 120000,
  });
  const before = await page.locator("#count").innerText();
  expect(Number(before)).toBeGreaterThan(0);
  const canvas = await page.locator("#viewport canvas").boundingBox();
  await page.mouse.click(
    canvas.x + canvas.width * 0.5,
    canvas.y + canvas.height * 0.5,
  );
  await expect(page.locator("#details dl")).toBeVisible();
  const skeletal = page.locator('#systems input[value="skeletal"]');
  await skeletal.uncheck();
  await expect(page.locator("#count")).not.toHaveText(before);
  await skeletal.check();
  await expect(page.locator("#count")).toHaveText(before);
  await page.locator("#structures button").first().click();
  await expect(page.locator("#details dl")).toBeVisible();
  await page.locator("#search").fill("this-structure-does-not-exist");
  await expect(page.locator("#count")).toHaveText("0");
  await expect(page.locator("#scene-status")).toContainText("No structures");
  await page.locator("#search").fill("");
  await page.locator("#opacity").fill("0.5");
  await page.locator("#opacity").dispatchEvent("input");
  await expect(page.locator("#opacity-value")).toHaveText("50%");
  await page.locator("#clip").fill("40");
  await page.locator("#clip").dispatchEvent("input");
  await expect(page.locator("#clip-value")).toHaveText("40%");
  await page.locator("#clip").fill("100");
  await page.locator("#clip").dispatchEvent("input");
  await page.locator("#posture").click();
  await expect(page.locator("#frame-label")).toContainText("Supine display");
  await page.locator("#posture").click();
  await page.locator("#tab-spectral").click();
  await expect(page.locator("#tab-spectral")).toHaveClass("active");
  await page.locator("#tab-phys").click();
  await page.locator("#opacity").fill("1");
  await page.locator("#opacity").dispatchEvent("input");
  await page.screenshot({ path: "test-results/workbench.png", fullPage: true });
  expect(errors).toEqual([]);
});
test("measured Laplace evidence and archived vascular playback use actual data", async ({
  page,
}) => {
  const errors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  await page.goto("/");
  await expect(page.locator("#variable")).toHaveValue(
    "ArterialPressure(mmHg)",
    { timeout: 30000 },
  );
  await page.locator("#tab-spectral").click();
  await page.locator("#spectral-run").selectOption("bidmc_01");
  await page.locator("#spectral-mode").selectOption("laplace");
  await page.locator("#sigma").selectOption("1");
  await expect(page.locator("#chart svg")).toBeVisible();
  await expect(page.locator("#chart-note")).toContainText("Finite-horizon");
  await page.locator("#model").selectOption("vascular-cerebral");
  await expect(page.locator("#play")).toBeEnabled({ timeout: 60000 });
  await page.locator("#play").click();
  await expect(page.locator("#time")).not.toHaveValue("0");
  await page.locator("#play").click();
  await page.locator("#flow-field").selectOption("pressure");
  await expect(page.locator("#flow-legend")).toContainText(
    "source units unconfirmed",
  );
  const legend = await page.locator("#flow-legend").innerText();
  await page.locator("#time").fill("150");
  await page.locator("#time").dispatchEvent("input");
  await expect(page.locator("#flow-legend")).toHaveText(legend);
  await page.screenshot({
    path: "test-results/vascular-laplace.png",
    fullPage: true,
  });
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(page.locator("#model")).toBeVisible();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBeTruthy();
  expect(errors).toEqual([]);
});
test("native scenario form submits bounded exercise and exposes the server run", async ({
  page,
}) => {
  test.skip(
    process.env.RUN_NATIVE_BROWSER !== "1",
    "Opt-in real native engine execution",
  );
  await page.goto("/");
  await expect(page.locator("#run-status")).not.toContainText("Checking");
  await page.locator("#scenario").selectOption("exercise");
  await page.locator("#duration").fill("2");
  const responsePromise = page.waitForResponse(
    (r) =>
      r.url().endsWith("/api/scenarios") && r.request().method() === "POST",
  );
  await page.locator("#run").click();
  const response = await responsePromise;
  expect(response.ok()).toBeTruthy();
  const run = await response.json();
  console.log("Native browser run:", run.id);
  await expect(page.locator("#run-status")).toContainText(run.id);
  await expect(page.locator("#run-status")).toContainText(
    /completed|complete|succeeded/,
    { timeout: 110000 },
  );
  await expect(page.locator("#chart svg")).toBeVisible();
});
test("OpenSim source family exposes attachment paths and unresolved wrapping", async ({
  page,
}) => {
  await page.goto("/");
  await page.locator("#model").selectOption("opensim-rajagopal");
  await expect(page.locator("#scene-status")).toHaveText("", {
    timeout: 60000,
  });
  await expect(page.locator("#frame-label")).toContainText("m display");
  await page.locator("#search").fill("addbrev_r");
  await page.locator("#structures button").first().click();
  await expect(page.locator("#details")).toContainText("attachment");
  await expect(page.locator("#details")).toContainText("Unsolved");
  await page.locator("#search").fill("");
  await page.screenshot({ path: "test-results/opensim.png", fullPage: true });
});
test("anatomy inspection remains available without a WebGL context", async ({
  page,
}) => {
  await page.addInitScript(() => {
    HTMLCanvasElement.prototype.getContext = () => null;
  });
  await page.goto("/");
  await expect(page.locator("#scene-status")).toContainText(
    "WebGL unavailable",
  );
  await expect(page.locator("#render-count")).toHaveText(
    "3D rendering unavailable",
  );
  await page.locator("#structures button").first().click();
  await expect(page.locator("#details dl")).toBeVisible();
});

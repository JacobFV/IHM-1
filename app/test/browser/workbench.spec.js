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
test("native scenario form submits selected-preset hemorrhage and saline with stop actions", async ({
  page,
}) => {
  test.skip(
    process.env.RUN_NATIVE_BROWSER !== "1",
    "Opt-in real native engine execution",
  );
  await page.goto("/");
  await expect(page.locator("#run-status")).not.toContainText("Checking");
  await page.locator("#patient").selectOption("StandardFemale");
  await expect(page.locator("#patient")).toHaveValue("StandardFemale");
  await page.locator("#patient").selectOption("StandardMale");
  await page.locator("#scenario").selectOption("hemorrhage_saline");
  await page.locator("#duration").fill("2");
  const responsePromise = page.waitForResponse(
    (r) =>
      r.url().endsWith("/api/scenarios") && r.request().method() === "POST",
  );
  await page.locator("#run").click();
  const response = await responsePromise;
  expect(response.ok()).toBeTruthy();
  const run = await response.json();
  expect(response.request().postDataJSON().patient).toBe("StandardMale");
  expect(
    response
      .request()
      .postDataJSON()
      .interventions.map((a) => a.kind),
  ).toEqual(["hemorrhage", "hemorrhage", "saline", "saline"]);
  console.log("Native browser run:", run.id);
  await expect(page.locator("#run-status")).toContainText(run.id);
  await expect(page.locator("#run-status")).toContainText(
    /completed|complete|succeeded/,
    { timeout: 110000 },
  );
  await expect(page.locator("#chart svg")).toBeVisible();
});
test("OpenSim source family exposes native wrapped paths and force provenance", async ({
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
  await expect(page.locator("#details")).toContainText("Native mechanics");
  await expect(page.locator("#details")).toContainText("Solved by source");
  await expect(page.locator("#details")).toContainText(
    "external force balance not solved",
  );
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
test("coverage, circuit balance and CFD uncertainty are exposed with measured and source-model context", async ({
  page,
}) => {
  const errors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  await page.goto("/");
  await expect(page.locator('#systems input[value="cardiac"]')).toBeChecked();
  await page.locator("summary").filter({ hasText: "System coverage" }).click();
  await expect(page.locator("#coverage-details")).toContainText(
    "21 system domains",
  );
  await page
    .locator("summary")
    .filter({ hasText: "Conservative exchange" })
    .click();
  await expect(page.locator("#coupling-details")).toContainText(
    "Maximum free-node flow residual",
  );
  await page
    .locator("summary")
    .filter({ hasText: "Vascular CFD audit" })
    .click();
  await expect(page.locator("#vascular-audit")).toContainText(
    "Simulation failed",
  );
  await expect(page.locator("#vascular-audit")).toContainText("24.9%");
  await page.locator("#trajectory-run").selectOption("reproductive");
  await page.locator("#variable").selectOption("E2");
  await expect(page.locator("#chart-note")).toContainText(
    "prescribed time input",
  );
  await expect(page.locator("#chart")).toContainText("day");
  await expect(page.locator("#chart svg")).toBeVisible();
  expect(errors).toEqual([]);
});

test("native BETSE cells expose voltage, ions, constant protein field and source-clock playback", async ({
  page,
}) => {
  const errors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  await page.goto("/");
  await page.locator("#model").selectOption("betse-tissue");
  await expect(page.locator("#play")).toBeEnabled();
  await expect(page.locator("#flow-legend")).toContainText("Vmem");
  await page.locator("#structures button").first().click();
  await expect(page.locator("#cell-values")).toContainText(
    "212 planar solver cells",
  );
  const canvas = page.locator("#viewport canvas"),
    box = await canvas.boundingBox();
  await canvas.click({ position: { x: box.width / 2, y: box.height / 2 } });
  await expect(page.locator("#cell-values strong")).toContainText("Cell");
  await page.locator("#flow-field").selectOption("proteins");
  await expect(page.locator("#flow-legend")).toContainText("135.0");
  await page.locator("#flow-field").selectOption("sodium");
  await expect(page.locator("#flow-legend")).toContainText("mol/m^3");
  await page.locator("#time").fill("33");
  await expect(page.locator("#time-value")).toContainText("0.034 s");
  await page.locator("#flow-field").selectOption("Vmem");
  await page.screenshot({ path: "test-results/betse-tissue.png" });
  expect(errors).toEqual([]);
});

test("CSF source trajectories retain coupling boundaries and environment overrides reach transport", async ({
  page,
}) => {
  await page.goto("/");
  for (const run of ["baseline", "native_map_driven", "hypotension"]) {
    await page.locator("#trajectory-run").selectOption(`csf:${run}`);
    await page.locator("#variable").selectOption("Pic_mmHg");
    await expect(page.locator("#chart svg")).toBeVisible();
    await expect(page.locator("#chart")).toContainText("mmHg");
    await expect(page.locator("#chart-note")).toContainText(
      run === "native_map_driven"
        ? "no ICP feedback"
        : "Separate literature model",
    );
  }
  await expect(page.locator("#ambient")).toHaveValue("");
  await expect(page.locator("#clothing")).toHaveValue("");
  let payload;
  await page.route("**/api/scenarios", async (route) => {
    if (route.request().method() !== "POST") return route.continue();
    payload = route.request().postDataJSON();
    await route.fulfill({
      status: 400,
      contentType: "application/json",
      body: JSON.stringify({ error: "Transport test: no native job executed" }),
    });
  });
  await page.locator("#ambient").fill("22");
  await page.locator("#clothing").fill("0");
  await page.locator("#run").click();
  await expect(page.locator("#run-status")).toContainText("Transport test");
  expect(payload.ambient_temperature_c).toBe(22);
  expect(payload.clothing_clo).toBe(0);
});

test("hour-scale spectra disclose sampling limits and exclude aliased pulse channels", async ({
  page,
}) => {
  await page.goto("/");
  await page.locator("#tab-spectral").click();
  await page.locator("#spectral-run").selectOption("native_hour_rest");
  await expect(page.locator("#chart-note")).toContainText("Nyquist 0.5 Hz");
  const names = await page.locator("#variable option").allTextContents();
  expect(names.some((n) => n === "ArterialPressure")).toBe(false);
  await expect(page.locator("#chart svg")).toBeVisible();
});

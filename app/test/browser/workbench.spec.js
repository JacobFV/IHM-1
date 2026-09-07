import { test, expect } from "@playwright/test";

test("the left column carries five sections and nothing asks for a source", async ({ page }) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto("/");
  await expect(page.locator("#viewport canvas")).toBeVisible();
  await expect(page.locator("#left-column .column-section > h2")).toHaveText([
    "Materialization", "Layers", "Clothing", "Environment", "Simulation",
  ]);
  // No source, candidate or execution pickers survive anywhere in the page.
  for (const id of ["#model", "#scene-owner", "#patient", "#engine-variant", "#scene-reconnect"])
    await expect(page.locator(id)).toHaveCount(0);
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  expect(errors).toEqual([]);
});

test("layers filter by member and expand to reveal members", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#layers .layer")).not.toHaveCount(0);
  // Members are listed only when the system is expanded.
  await expect(page.locator("#layers .members")).toHaveCount(0);
  await page.locator("#layer-search").fill("lung");
  await expect(page.locator("#layers .layer-row > label > span")).toHaveText(["respiratory"]);
  await expect(page.locator("#layers .members .check-row")).not.toHaveCount(0);
  await page.locator("#layer-search").fill("");
  await expect(page.locator("#layers .members")).toHaveCount(0);
  await page.locator("#layers .disclose").first().click();
  await expect(page.locator("#layers .members")).toHaveCount(1);
});

test("clothing and environment appear only for the whole body", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#clothing-section")).toBeVisible();
  await expect(page.locator("#environment-section")).toBeVisible();
  await page.locator("#materialization option[value^='conforming-domain-']").first().waitFor({ state: "attached" });
  const domain = await page.locator("#materialization option[value^='conforming-domain-']").first().getAttribute("value");
  await page.locator("#materialization").selectOption(domain);
  await expect(page.locator("#clothing-section")).toBeHidden();
  await expect(page.locator("#environment-section")).toBeHidden();
  await expect(page.locator("[data-pane='domain'] .check-row")).not.toHaveCount(0, { timeout: 240000 });
  await page.locator("#materialization").selectOption("body");
  await expect(page.locator("#clothing-section")).toBeVisible();
});

test("the viewport carries only the gimbal, the transport and two column icons", async ({ page }) => {
  await page.goto("/");
  // One human figure carrying three selectable planes, not a row of buttons.
  await expect(page.locator("#gimbal .figure")).toHaveCount(1);
  await expect(page.locator("#gimbal [data-plane]")).toHaveCount(3);
  await expect(page.locator("#speed")).toHaveText("1x");
  await page.locator("#speed").click();
  await expect(page.locator("#speed-menu button")).toHaveText(["0.25x", "0.5x", "1x", "2x", "5x", "10x"]);
  await page.locator("#speed-menu button").getByText("5x", { exact: true }).click();
  await expect(page.locator("#speed")).toHaveText("5x");
  // The composed body-state legend is gone, not shortened.
  await expect(page.locator("#flow-legend")).toHaveCount(0);
  await expect(page.getByText("tissue transforms")).toHaveCount(0);
  await page.locator("#toggle-left").click();
  await expect(page.locator("#left-column")).toBeHidden();
  await page.locator("#toggle-right").click();
  await expect(page.locator("#pane-column")).toBeHidden();
  await expect(page.locator("#gimbal")).toBeVisible();
  await expect(page.locator("#transport")).toBeVisible();
  await page.locator("#toggle-left").click();
  await expect(page.locator("#left-column")).toBeVisible();
});

test("monitors can be removed and chosen again", async ({ page }) => {
  await page.goto("/");
  await page.locator("#add-pane").waitFor({ timeout: 60000 });
  await expect(page.locator("[data-pane='selection']")).toBeVisible();
  await page.locator("[data-pane='selection'] .pane-dismiss").click();
  await expect(page.locator("[data-pane='selection']")).toBeHidden();
  await page.locator("#add-pane").click();
  const entry = page.locator("#pane-picker button").getByText("Selection", { exact: true });
  await expect(entry).toHaveCount(1);
  await entry.click();
  await expect(page.locator("[data-pane='selection']")).toBeVisible();
});

test("clicking a structure identifies it and answers where it came from", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  const canvas = await page.locator("#viewport canvas").boundingBox();
  await page.mouse.click(canvas.x + canvas.width * 0.5, canvas.y + canvas.height * 0.42);
  await expect(page.locator("[data-pane='selection'] h3")).toBeVisible();
  await expect(page.locator("#details .structure-provenance")).toContainText("Where this came from", { timeout: 60000 });
  await expect(page.locator("#details .structure-provenance")).toContainText("Originating dataset");
});

import { test, expect } from "@playwright/test";

test("the left column carries five sections and nothing asks for a source", async ({ page }) => {
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto("/");
  await expect(page.locator("#scene")).toBeVisible();
  await expect(page.locator("#left-column .column-section h2")).toHaveText([
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

test("the environment section tells tiles, configuration and objects apart", async ({ page }) => {
  await page.goto("/");
  const tile = (id) => page.locator(`#environment-tiles [data-tile="${id}"]`);
  const config = (slot) => page.locator(`#environment-tiles select[data-slot="${slot}"]`);
  await expect(tile("bed")).toHaveAttribute("aria-pressed", "true", { timeout: 60000 });
  // Environments and scenes are tiles; components never are.
  await expect(tile("bedroom")).toBeVisible();
  await expect(tile("mattress-soft")).toHaveCount(0);
  // Configuration is a dropdown that appears only when its requirements hold.
  await expect(config("ambient_thermal")).toBeVisible();
  await expect(config("bed_support_model")).toBeVisible();
  await expect(config("mattress_material")).toHaveCount(0);
  await config("bed_support_model").selectOption("bed-support-skin-quadrature");
  await expect(config("mattress_material")).toBeVisible();
  await config("mattress_material").selectOption("mattress-soft");
  await expect(config("mattress_material")).toHaveValue("mattress-soft");
  // Leaving the bed withdraws everything that depended on it.
  await tile("studio").click();
  await expect(config("mattress_material")).toHaveCount(0);
  await expect(config("bed_support_model")).toHaveCount(0);
  await expect(tile("bedroom")).toBeDisabled();
  // Ambient air is independent of the mechanical environment.
  await expect(config("ambient_thermal")).toBeVisible();
});

test("scene objects are additive inserts, not toggles", async ({ page }) => {
  await page.goto("/");
  const insert = (id) => page.locator(`#environment-tiles [data-insert="${id}"]`);
  await expect(insert("ball-small")).toBeVisible({ timeout: 60000 });
  await expect(page.locator("#environment-tiles .instance-list li")).toHaveCount(0);
  // The same object can be inserted more than once.
  await insert("ball-small").click();
  await insert("ball-small").click();
  await insert("pillow").click();
  await expect(page.locator("#environment-tiles .instance-list li")).toHaveCount(3);
  await page.locator("#environment-tiles .instance-list li .instance-remove").first().click();
  await expect(page.locator("#environment-tiles .instance-list li")).toHaveCount(2);
  // An instance whose environment requirement lapses leaves with it; the ball
  // requires nothing, so it stays.
  await page.locator('#environment-tiles [data-tile="studio"]').click();
  await expect(page.locator("#environment-tiles .instance-list li")).toHaveCount(1);
});

test("the viewport carries only the gimbal, the transport and two column icons", async ({ page }) => {
  await page.goto("/");
  // A real 3D widget carrying this body, not a drawing: its own canvas, and it
  // turns with the main camera.
  await expect(page.locator("canvas#gimbal")).toBeVisible();
  await expect(page.locator("#gimbal")).toHaveAttribute("data-selected", "coronal");
  const before = await page.locator("#gimbal").screenshot();
  const box = await page.locator("#scene").boundingBox();
  await page.mouse.move(box.x + box.width * 0.5, box.y + box.height * 0.5);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.5 + 220, box.y + box.height * 0.5 + 40, { steps: 10 });
  await page.mouse.up();
  await page.waitForTimeout(800);
  expect((await page.locator("#gimbal").screenshot()).equals(before)).toBe(false);
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

test("all three planes are reachable on the gimbal and snap the camera", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("canvas#gimbal")).toBeVisible();
  const box = await page.locator("#scene").boundingBox();
  const rect = await page.locator("#gimbal").boundingBox();
  // From an axis-aligned view two planes are exactly edge-on, which is correct
  // for a view gizmo, so orbit to a general orientation first.
  await page.mouse.move(box.x + box.width * 0.5, box.y + box.height * 0.5);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.5 + 190, box.y + box.height * 0.5 + 65, { steps: 10 });
  await page.mouse.up();
  await page.waitForTimeout(1500);
  // Every plane is found by real hit testing, which the widget reports as the
  // pointer crosses it; nothing here reaches into the scene graph.
  // The scan runs in one round trip: each point is a real pointermove on the
  // canvas, so the widget's own hit testing decides, exactly as for a cursor.
  const found = await page.evaluate(() => {
    const canvas = document.getElementById("gimbal");
    const box = canvas.getBoundingClientRect();
    const hits = {};
    for (let i = 1; i < 16; i++)
      for (let j = 1; j < 16; j++) {
        const clientX = box.left + (box.width * i) / 16, clientY = box.top + (box.height * j) / 16;
        canvas.dispatchEvent(new PointerEvent("pointermove", { clientX, clientY, bubbles: true }));
        const hover = canvas.dataset.hover;
        if (hover && !hits[hover]) hits[hover] = { x: clientX, y: clientY };
      }
    return hits;
  });
  expect(Object.keys(found).sort()).toEqual(["coronal", "sagittal", "transverse"]);
  // Clicking one snaps the main camera to that view. Only one snap is asserted:
  // a snap leaves the other two edge-on, and re-probing races camera damping.
  const plane = "transverse";
  await page.mouse.move(found[plane].x, found[plane].y);
  await expect(page.locator("#gimbal")).toHaveAttribute("data-hover", plane);
  const before = await page.locator("#scene").screenshot();
  await page.mouse.click(found[plane].x, found[plane].y);
  await expect(page.locator("#gimbal")).toHaveAttribute("data-selected", plane);
  await page.waitForTimeout(1500);
  expect((await page.locator("#scene").screenshot()).equals(before)).toBe(false);
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

test("every left section collapses to its title and stays collapsed", async ({ page }) => {
  await page.goto("/");
  const sections = page.locator("#left-column .column-section");
  await expect(sections).toHaveCount(5);
  await expect(sections.locator(".section-head .disclose")).toHaveCount(5);
  const layers = page.locator("[data-section='layers']");
  const open = (await layers.boundingBox()).height;
  await layers.locator(".section-head .disclose").click();
  await expect(layers).toHaveAttribute("data-collapsed", "true");
  await expect(layers.locator(".section-body")).toBeHidden();
  expect((await layers.boundingBox()).height).toBeLessThan(open);
  await page.reload();
  await expect(page.locator("[data-section='layers']")).toHaveAttribute("data-collapsed", "true");
  await page.locator("[data-section='layers'] .section-head .disclose").click();
  await expect(page.locator("[data-section='layers']")).toHaveAttribute("data-collapsed", "false");
});

test("clothing and environment tiles share one compact footprint", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#environment-tiles .tile").first()).toBeVisible({ timeout: 60000 });
  const box = (locator) => locator.first().boundingBox();
  const garment = await box(page.locator("#clothing-tiles .tile"));
  const environment = await box(page.locator("#environment-tiles .tile"));
  expect(Math.abs(garment.width - environment.width)).toBeLessThan(2);
  expect(Math.abs(garment.height - environment.height)).toBeLessThan(2);
  // Two per row in both grids.
  const rows = await page.locator("#environment-tiles .tile-row").first().evaluate((el) =>
    getComputedStyle(el).gridTemplateColumns.split(" ").length);
  expect(rows).toBe(2);
  // No section scrolls inside itself; the column scrolls as a whole.
  const inner = await page.locator("#left-column .section-body, #layers, .members, #environment-tiles")
    .evaluateAll((nodes) => nodes.filter((n) => {
      const style = getComputedStyle(n);
      return ["auto", "scroll"].includes(style.overflowY) || style.maxHeight !== "none";
    }).length);
  expect(inner).toBe(0);
  // The picture is not desaturated, tinted or overlaid anywhere in our styling.
  const picture = await page.locator("#environment-tiles .tile img").first().evaluate((el) => {
    const style = getComputedStyle(el);
    const parent = getComputedStyle(el.parentElement);
    return { filter: style.filter, opacity: style.opacity, mix: style.mixBlendMode, parentFilter: parent.filter };
  });
  expect(picture.filter).toBe("none");
  expect(picture.parentFilter).toBe("none");
  expect(picture.opacity).toBe("1");
  expect(picture.mix).toBe("normal");
});

test("every pane collapses to its title and stays collapsed across a reload", async ({ page }) => {
  await page.goto("/");
  await page.locator("#add-pane").waitFor({ timeout: 60000 });
  const panes = page.locator("#pane-column .pane");
  // Uniformly: every pane has a disclosure and a remove, both always available.
  await expect(panes.locator(".disclose")).toHaveCount(await panes.count());
  await expect(panes.locator(".pane-dismiss")).toHaveCount(await panes.count());
  const pane = page.locator("[data-pane='selection']");
  const open = (await pane.boundingBox()).height;
  await pane.locator(".disclose").click();
  await expect(pane).toHaveAttribute("data-collapsed", "true");
  await expect(pane.locator(".pane-body")).toBeHidden();
  await expect(pane.locator(".pane-dismiss")).toBeVisible();
  const shut = (await pane.boundingBox()).height;
  expect(shut).toBeLessThan(open);
  await page.reload();
  await expect(page.locator("[data-pane='selection']")).toHaveAttribute("data-collapsed", "true");
  await page.locator("[data-pane='selection'] .disclose").click();
  await expect(page.locator("[data-pane='selection']")).toHaveAttribute("data-collapsed", "false");
});

test("clicking a structure identifies it and answers where it came from", async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  const canvas = await page.locator("#scene").boundingBox();
  await page.mouse.click(canvas.x + canvas.width * 0.5, canvas.y + canvas.height * 0.42);
  await expect(page.locator("[data-pane='selection'] h3")).toBeVisible();
  await expect(page.locator("#details .structure-provenance")).toContainText("Where this came from", { timeout: 60000 });
  await expect(page.locator("#details .structure-provenance")).toContainText("Originating dataset");
});

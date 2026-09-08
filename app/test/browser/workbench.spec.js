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
  // Ambient air is independent of the mechanical environment.
  await expect(config("ambient_thermal")).toBeVisible();
  // A scene that names one environment is still offered from another, and says
  // what pressing it will also do, because a requirement is a route rather than
  // a wall. Configuration, which offers a real choice, stays withdrawn.
  await expect(tile("bedroom")).toBeEnabled();
  await expect(tile("bedroom")).toHaveAttribute("title", "Also selects Bed");
  await tile("bedroom").click();
  await expect(tile("bed")).toHaveAttribute("aria-pressed", "true");
  await expect(tile("bedroom")).toHaveAttribute("aria-pressed", "true");
  await expect(config("bed_support_model")).toBeVisible();
});

test("scene objects are additive inserts, not toggles", async ({ page }) => {
  await page.goto("/");
  // One plus opens a menu of the insertable objects; the buttons live inside it.
  const plus = page.locator("#insert-object");
  const insert = async (id) => {
    await plus.click();
    await page.locator(`#environment-tiles [data-insert="${id}"]`).click();
  };
  await expect(plus).toBeVisible({ timeout: 60000 });
  await expect(page.locator("#insert-menu")).toBeHidden();
  await expect(page.locator("#environment-tiles .instance-list li")).toHaveCount(0);
  // The same object can be inserted more than once.
  await insert("ball-small");
  await insert("ball-small");
  await insert("pillow");
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
  // All three are reachable, asserted above. For the snap, take whichever plane
  // the widget reports under the pointer and which is not already selected:
  // which plane sits at a given pixel depends on an orbit that is not
  // bit-deterministic, so naming one in advance would test the drag, not the
  // gimbal. Scan and confirm in one loop, because the widget keeps turning.
  const already = await page.locator("#gimbal").getAttribute("data-selected");
  const before = await page.locator("#scene").screenshot();
  let plane = null;
  for (let attempt = 0; attempt < 12 && !plane; attempt++) {
    const at = await page.evaluate((skip) => {
      const canvas = document.getElementById("gimbal");
      const box = canvas.getBoundingClientRect();
      for (let i = 1; i < 24; i++)
        for (let j = 1; j < 24; j++) {
          const clientX = box.left + (box.width * i) / 24, clientY = box.top + (box.height * j) / 24;
          canvas.dispatchEvent(new PointerEvent("pointermove", { clientX, clientY, bubbles: true }));
          const hit = canvas.dataset.hover;
          if (hit && hit !== skip) return { x: clientX, y: clientY, hit };
        }
      return null;
    }, already);
    if (!at) { await page.waitForTimeout(250); continue; }
    await page.mouse.move(at.x, at.y);
    if ((await page.locator("#gimbal").getAttribute("data-hover")) !== at.hit) continue;
    // Press where the pointer already is, in the same settled state the hover
    // was confirmed against; the widget turns with the camera between probes.
    await page.mouse.down();
    await page.mouse.up();
    if ((await page.locator("#gimbal").getAttribute("data-selected")) === at.hit) plane = at.hit;
  }
  expect(["coronal", "sagittal", "transverse"]).toContain(plane);
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

test("an equal drag turns the camera equally in every direction, and the arc over the head does not stop", async ({ page }) => {
  // Every reading is a real drag on the canvas and a settled camera read back,
  // and the software rasteriser makes both slow; the exhaustive sweep across
  // five orientations lives in app/test/orbit-isotropy.mjs.
  test.setTimeout(300000);
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  const box = await page.locator("#scene").boundingBox();
  const centre = { x: box.x + box.width * 0.5, y: box.y + box.height * 0.5 };
  // A drag of one canvas height is one full turn, so the expected angle is a
  // number this test states rather than a value read off the thing measured.
  const canvasHeight = await page.locator("#scene").evaluate((node) => node.clientHeight);
  const PIXELS = 140;
  const expected = (360 * PIXELS) / canvasHeight;
  const angle = (a, b) => {
    const d = Math.min(1, Math.max(-1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]));
    return (Math.acos(d) * 180) / Math.PI;
  };
  // Damping spreads a drag over many frames; this runs the same update the
  // render loop runs until it is spent, so a reading is of a finished camera.
  const settle = () => page.evaluate(() => globalThis.__ihmCamera.settle(900));
  const drag = async (dx, dy) => {
    await page.mouse.move(centre.x, centre.y);
    await page.mouse.down();
    await page.mouse.move(centre.x + dx, centre.y + dy, { steps: 4 });
    await page.mouse.up();
    return settle();
  };

  // Three starts, two drag directions each: the same pixels must buy the same
  // degrees whichever way they go and wherever they start. The transverse snap
  // is the telling one -- it is exactly where a turntable's azimuth degenerates,
  // and a drag there used to move the camera by a ten-thousandth of a degree.
  const readings = [];
  for (const plane of ["coronal", "sagittal", "transverse"])
    for (const [dx, dy] of [[PIXELS, 0], [0, PIXELS]]) {
      await page.evaluate((p) => globalThis.__ihmCamera.snap(p), plane);
      const before = await settle();
      const after = await drag(dx, dy);
      readings.push({ plane, dx, dy, deg: angle(before.direction, after.direction) });
    }
  for (const reading of readings)
    expect.soft(reading.deg, `${reading.plane} drag ${reading.dx},${reading.dy}`).toBeCloseTo(expected, 1);
  const degrees = readings.map((r) => r.deg);
  expect(Math.max(...degrees) / Math.min(...degrees)).toBeLessThan(1.02);

  // The owner's complaint, as an assertion: hold a vertical drag from the front
  // view and the camera travels the coronal plane, over the head, past the back
  // and round. A turntable stops on the pole after a quarter turn.
  await page.evaluate(() => globalThis.__ihmCamera.snap("coronal"));
  let previous = await settle();
  const first = previous;
  let travelled = 0, overhead = false, behind = false;
  // Four pulls of 200 px is 288 degrees: far enough round to pass behind the
  // body and over the head, and deliberately short of a full turn so the snap
  // below has somewhere to come back from.
  for (let pull = 0; pull < 4; pull++) {
    const now = await drag(0, -200);
    travelled += angle(previous.direction, now.direction);
    previous = now;
    if (now.direction[1] > 0.8) overhead = true;
    if (now.direction[2] < -0.8) behind = true;
  }
  expect(travelled).toBeCloseTo((360 * 4 * 200) / canvasHeight, 0);
  expect(overhead, "the camera never reached the view over the head").toBe(true);
  expect(behind, "the camera never carried on past the head to the back").toBe(true);

  // A snap still establishes an upright view from wherever the sweep ended, and
  // the gimbal still tracks the camera it mirrors.
  // The sweep left the camera 72 degrees short of home, so the gimbal, which
  // mirrors the camera every frame, must look different after the snap.
  const gimbalBefore = await page.locator("#gimbal").screenshot();
  await page.evaluate(() => globalThis.__ihmCamera.snap("coronal"));
  const home = await settle();
  expect(angle(home.direction, [0, 0, 1])).toBeLessThan(0.01);
  expect(angle(home.up, [0, 1, 0])).toBeLessThan(0.01);
  expect(angle(first.direction, home.direction)).toBeLessThan(0.01);
  await page.waitForTimeout(800);
  expect((await page.locator("#gimbal").screenshot()).equals(gimbalBefore)).toBe(false);
  await expect(page.locator("#gimbal")).toHaveAttribute("data-selected", "coronal");
  expect(errors).toEqual([]);
});

test("every system row has an opacity control that changes only that system", async ({ page }) => {
  test.setTimeout(300000);
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  const rows = page.locator("#layers .layer-row");
  const count = await rows.count();
  expect(count).toBeGreaterThan(3);
  // One control on every system row, and the row at rest is otherwise what it
  // was: a triangle, a checkbox, a swatch, a name and a count. No inline slider.
  await expect(page.locator("#layers .layer-row .layer-settings")).toHaveCount(count);
  await expect(page.locator("#layers .layer-row input[type=range]")).toHaveCount(0);
  await expect(page.locator("#layer-opacity")).toBeHidden();

  const painted = (system) => page.evaluate((s) => globalThis.__ihmLayers.opacity(s), system);
  const control = (system) => page.locator(`#layers .layer-settings[data-system="${system}"]`);
  const systems = await page.locator("#layers .layer-settings").evaluateAll((n) => n.map((x) => x.dataset.system));
  expect(systems).toContain("integumentary");
  // A system whose geometry is actually in the scene, so a change to it is a
  // change to the body rather than to a number nothing paints.
  let other = null;
  for (const system of systems)
    if (system !== "integumentary" && (await painted(system)) !== null) { other = system; break; }
  expect(other, "no other system is currently painted").toBeTruthy();

  // Opening one and then another leaves exactly one popover open.
  await control("integumentary").click();
  await expect(page.locator("#layer-opacity")).toBeVisible();
  await expect(page.locator("#layer-opacity")).toHaveAttribute("data-system", "integumentary");
  await expect(control("integumentary")).toHaveAttribute("aria-expanded", "true");
  await control(other).click();
  await expect(page.locator("#layer-opacity")).toHaveAttribute("data-system", other);
  await expect(control("integumentary")).toHaveAttribute("aria-expanded", "false");
  await expect(page.locator("#layer-opacity")).toHaveCount(1);

  // The slider moves the body, not just the readout, and only its own system.
  const untouched = await painted("integumentary");
  await page.locator("#layer-opacity-range").fill("35");
  await expect(page.locator("#layer-opacity-value")).toHaveText("35%");
  expect(await painted(other)).toBeCloseTo(0.35, 2);
  expect(await painted("integumentary")).toBeCloseTo(untouched, 5);
  // Opacity is not visibility: the checkbox is untouched.
  await expect(page.locator(`#layers .layer-row input[type=checkbox][value="${other}"]`)).toBeChecked();

  // Escape closes it; so does a click outside.
  await page.keyboard.press("Escape");
  await expect(page.locator("#layer-opacity")).toBeHidden();
  await control(other).click();
  await expect(page.locator("#layer-opacity")).toBeVisible();
  await page.locator("#layer-search").click();
  await expect(page.locator("#layer-opacity")).toBeHidden();

  // A row near the bottom of a tall column opens its popover on screen.
  const last = page.locator("#layers .layer-settings").last();
  await last.scrollIntoViewIfNeeded();
  await last.click();
  const placed = await page.locator("#layer-opacity").boundingBox();
  const viewportHeight = await page.evaluate(() => window.innerHeight);
  expect(placed.y).toBeGreaterThanOrEqual(0);
  expect(placed.y + placed.height).toBeLessThanOrEqual(viewportHeight);
  await page.keyboard.press("Escape");

  // The value survives a reload, and still does not touch visibility.
  await page.reload();
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  expect(await painted(other)).toBeCloseTo(0.35, 2);
  await control(other).click();
  await expect(page.locator("#layer-opacity-value")).toHaveText("35%");
  await expect(page.locator(`#layers .layer-row input[type=checkbox][value="${other}"]`)).toBeChecked();
});

test("the Layers section picks a tissue colour palette, and palette and skin tone are separate axes", async ({ page }) => {
  test.setTimeout(300000);
  const errors = [];
  page.on("pageerror", (error) => errors.push(error.message));
  await page.goto("/");
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });

  // The control belongs to Layers and is the same gear-and-popover idiom as the
  // per-system opacity control, not a second kind of thing.
  await expect(page.locator("[data-section='layers'] #palette-row #palette-settings")).toBeVisible();
  await expect(page.locator("#tissue-palette")).toBeHidden();
  // The build default is what the app opens at, and the popover says so in the
  // server's own words rather than in a claim written here.
  await expect(page.locator("#palette-current")).toHaveText("Didactic");

  // Three structures whose colour the reader can see: the whole-body skin mesh,
  // the lip, and the gingiva. Didactic paints the last two the same because they
  // are one system; a tissue palette must not.
  const painted = (id) => page.evaluate((i) => globalThis.__ihmLayers.colour(i), id);
  const skin = "body-bp3d-FJ2810", lip = "body-bp3d-FJ2814", gingiva = "body-bp3d-FJ1252";
  const skinDidactic = await painted(skin);
  expect(skinDidactic).toMatch(/^#[0-9a-f]{6}$/);
  expect(await painted(lip)).toBe(await painted(gingiva));

  await page.locator("#palette-settings").click();
  await expect(page.locator("#tissue-palette")).toBeVisible();
  await expect(page.locator("#palette-settings")).toHaveAttribute("aria-expanded", "true");
  await expect(page.locator("#palette-choice option")).toHaveText([
    "Didactic", "Realistic (fresh, in vivo)", "Evidence tier",
  ]);
  // Skin tone is published beside the palettes, not inside one, and is its own
  // control. Didactic declares no skin tone, so the control says so and is off.
  await expect(page.locator("#skin-tone-choice")).toBeDisabled();
  await expect(page.locator("#skin-tone-note")).toContainText("declares no skin tone");
  await expect(page.locator("#palette-policy")).toContainText("didactic stays the default");

  // The screenshot is taken before the click that is supposed to change it.
  const beforePalette = await page.locator("#scene").screenshot();
  await page.locator("#palette-choice").selectOption("realistic");
  await expect(page.locator("#palette-current")).toHaveText("Realistic (fresh, in vivo)");
  await expect(page.locator("#palette-note")).toContainText("not a formalin-fixed cadaver");
  // A recolour, not a reload: the body is repainted where it stands.
  // Both sides are re-read on every poll; comparing against one value read once
  // here would compare the body against whatever it was before the click.
  await expect.poll(async () => (await painted(lip)) === (await painted(gingiva)),
    { timeout: 30000 }).toBe(false);
  await page.waitForTimeout(1200);
  expect((await page.locator("#scene").screenshot()).equals(beforePalette)).toBe(false);

  // The skin tone moves the skin the palette declares a tone for, and leaves
  // every tissue whose colour was measured on its own exactly where it was.
  await expect(page.locator("#skin-tone-choice")).toBeEnabled();
  const skinNeutral = await painted(skin), lipRealistic = await painted(lip);
  const beforeTone = await page.locator("#scene").screenshot();
  await page.locator("#skin-tone-choice").selectOption("ita_dark");
  await expect.poll(() => painted(skin), { timeout: 30000 }).not.toBe(skinNeutral);
  const skinDark = await painted(skin);
  expect(await painted(lip)).toBe(lipRealistic);
  await page.waitForTimeout(1200);
  expect((await page.locator("#scene").screenshot()).equals(beforeTone)).toBe(false);

  // Choosing a palette that declares no skin tone turns the tone control off
  // again without losing the tone that was chosen.
  await page.locator("#palette-choice").selectOption("evidence");
  await expect(page.locator("#palette-current")).toHaveText("Evidence tier");
  await expect(page.locator("#skin-tone-choice")).toBeDisabled();
  await expect(page.locator("#skin-tone-choice")).toHaveValue("ita_dark");
  await page.locator("#palette-choice").selectOption("realistic");
  await expect.poll(() => painted(skin), { timeout: 30000 }).toBe(skinDark);

  // Escape closes it, a click outside closes it, and opening it closes an open
  // opacity popover: one popover at a time, as before.
  await page.keyboard.press("Escape");
  await expect(page.locator("#tissue-palette")).toBeHidden();
  await page.locator("#layers .layer-settings").first().click();
  await expect(page.locator("#layer-opacity")).toBeVisible();
  await page.locator("#palette-settings").click();
  await expect(page.locator("#layer-opacity")).toBeHidden();
  await expect(page.locator("#tissue-palette")).toBeVisible();
  await page.locator("#layer-search").click();
  await expect(page.locator("#tissue-palette")).toBeHidden();

  // Both axes survive a reload, and no geometry had to be refetched to apply
  // them: the reloaded body comes up already painted.
  await page.reload();
  await expect(page.locator("#scene-status")).toHaveText("", { timeout: 240000 });
  await expect(page.locator("#palette-current")).toHaveText("Realistic (fresh, in vivo)");
  await expect.poll(() => painted(skin), { timeout: 60000 }).toBe(skinDark);
  expect(skinDark).not.toBe(skinDidactic);
  await page.locator("#palette-settings").click();
  await expect(page.locator("#palette-choice")).toHaveValue("realistic");
  await expect(page.locator("#skin-tone-choice")).toHaveValue("ita_dark");
  expect(errors).toEqual([]);
});

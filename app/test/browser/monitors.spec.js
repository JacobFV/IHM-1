import { test, expect } from "@playwright/test";

test("full-height scene with configurable monitor cards and preserved widgets", async ({
  page,
}) => {
  await page.goto("/");
  await expect(page.locator("#viewport canvas")).toBeVisible();
  await expect(page.locator("body > #app > header")).toHaveCount(0);
  await expect(page.locator("#monitor-stack")).toBeVisible();
  expect(
    (await page.locator("#viewport").boundingBox()).height,
  ).toBeGreaterThan(950);
  await page.locator("#duration").fill("17");
  await page
    .getByRole("button", { name: "Remove Experiments", exact: true })
    .click();
  await expect(page.locator("#scenario-form")).toBeHidden();
  await page.getByRole("button", { name: "New pane", exact: true }).click();
  await page
    .getByRole("searchbox", { name: "Search monitor catalog" })
    .fill("protocol");
  await page
    .getByRole("button", { name: "Add Experiments", exact: true })
    .click();
  await expect(page.locator("#duration")).toHaveValue("17");
  await page
    .getByRole("button", { name: "Collapse Signal graph", exact: true })
    .click();
  await expect(page.locator("#chart")).toBeHidden();
  await page
    .getByRole("button", { name: "Expand Signal graph", exact: true })
    .click();
  await expect(page.locator("#chart")).toBeVisible();
  await page
    .getByRole("button", { name: "Configure Signal graph", exact: true })
    .click();
  await page.getByRole("checkbox", { name: "Compact Signal graph" }).check();
  await expect(page.locator('[data-monitor="signals"]')).toHaveAttribute(
    "data-density",
    "compact",
  );
  expect(
    await page.evaluate(
      () =>
        document.documentElement.scrollHeight <= innerHeight &&
        document.documentElement.scrollWidth <= innerWidth,
    ),
  ).toBe(true);
  await page.reload();
  await expect(page.locator('[data-monitor="signals"]')).toHaveAttribute(
    "data-density",
    "compact",
  );
});

test("catalog uses overlapping tags and keyboard dismissal", async ({
  page,
}) => {
  await page.goto("/");
  await page.getByRole("button", { name: "New pane", exact: true }).click();
  await page
    .getByRole("searchbox", { name: "Search monitor catalog" })
    .fill("mechanics");
  await expect(
    page.getByRole("button", { name: "Add Local contact", exact: true }),
  ).toBeVisible();
  await page.keyboard.press("Escape");
  await expect(
    page.getByRole("dialog", { name: "Monitor catalog" }),
  ).toBeHidden();
  await expect(
    page.getByRole("button", { name: "New pane", exact: true }),
  ).toBeFocused();
});

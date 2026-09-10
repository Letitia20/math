import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const [templatePath, payloadPath, outputPath] = process.argv.slice(2);
if (!templatePath || !payloadPath || !outputPath) {
  throw new Error("Usage: node tools/build_result2.mjs TEMPLATE PAYLOAD OUTPUT");
}

const payload = JSON.parse(await fs.readFile(payloadPath, "utf8"));
const rowCount = payload.time_s.length;
const radiusCount = payload.radius_cm.length;
if (rowCount !== 10800 || radiusCount !== 21) {
  throw new Error(`Expected 10800 times and 21 radii; got ${rowCount} and ${radiusCount}`);
}

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(templatePath));
const sheetData = [
  ["温度", payload.temperature_c],
  ["水分浓度", payload.moisture_concentration],
];

for (const [sheetName, field] of sheetData) {
  if (field.length !== rowCount || field.some((row) => row.length !== radiusCount)) {
    throw new Error(`${sheetName} has an invalid result matrix shape`);
  }
  const sheet = workbook.worksheets.getItem(sheetName);
  const matrix = [
    ["时间\\到药材中心的距离", ...payload.radius_cm],
    ...payload.time_s.map((time, index) => [time, ...field[index]]),
  ];
  const outputRange = sheet.getRangeByIndexes(0, 0, rowCount + 1, radiusCount + 1);
  outputRange.values = matrix;
  outputRange.format.font = { name: "宋体", size: 10 };
  outputRange.format.horizontalAlignment = "center";
  outputRange.format.verticalAlignment = "center";
  sheet.getRangeByIndexes(1, 0, rowCount, 1).format.numberFormat = "0";
  sheet.getRangeByIndexes(1, 1, rowCount, radiusCount).format.numberFormat = "0.0000";
  sheet.getRange("A:A").format.columnWidth = 24;
  sheet.getRange("B:V").format.columnWidth = 11;
  sheet.getRange("1:1").format.rowHeight = 22;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
}

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(outputPath);

for (const [sheetName] of sheetData) {
  for (const range of ["A1:V6", "A10798:V10801"]) {
    const check = await workbook.inspect({
      kind: "region",
      sheetId: sheetName,
      range,
      maxChars: 7000,
      tableMaxRows: 8,
      tableMaxCols: 22,
    });
    console.log(check.ndjson);
  }
}

await fs.mkdir("outputs/previews", { recursive: true });
for (const [sheetName] of sheetData) {
  const preview = await workbook.render({
    sheetName,
    range: "A1:V26",
    scale: 1.25,
    format: "png",
  });
  await fs.writeFile(
    `outputs/previews/problem2-${sheetName}.png`,
    new Uint8Array(await preview.arrayBuffer()),
  );
}

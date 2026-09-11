import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const [templatePath, payloadPath, outputPath] = process.argv.slice(2);
if (!templatePath || !payloadPath || !outputPath) {
  throw new Error("Usage: node tools/build_result3.mjs TEMPLATE PAYLOAD OUTPUT");
}

const payload = JSON.parse(await fs.readFile(payloadPath, "utf8"));
const rowCount = payload.time_s.length;
const radiusCount = payload.radius_cm.length;
if (rowCount < 2 || radiusCount !== 21) {
  throw new Error(`Expected at least 2 times and 21 radii; got ${rowCount} and ${radiusCount}`);
}
if (
  payload.moisture_concentration.length !== rowCount ||
  payload.moisture_concentration.some((row) => row.length !== radiusCount)
) {
  throw new Error("Moisture result matrix has an invalid shape");
}

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(templatePath));
const sheet = workbook.worksheets.getItem("Sheet1");
const matrix = [
  ["时间\\到药材中心的距离", ...payload.radius_cm],
  ...payload.time_s.map((time, index) => [
    time,
    ...payload.moisture_concentration[index],
  ]),
];
const outputRange = sheet.getRangeByIndexes(0, 0, rowCount + 1, radiusCount + 1);
outputRange.values = matrix;
outputRange.format.font = { name: "宋体", size: 10 };
outputRange.format.horizontalAlignment = "center";
outputRange.format.verticalAlignment = "center";
sheet.getRangeByIndexes(1, 0, rowCount, 1).format.numberFormat = "0.######";
sheet.getRangeByIndexes(1, 1, rowCount, radiusCount).format.numberFormat = "0.0000";
sheet.getRange("A:A").format.columnWidth = 24;
sheet.getRange("B:V").format.columnWidth = 11;
sheet.getRange("1:1").format.rowHeight = 22;
sheet.freezePanes.freezeRows(1);
sheet.freezePanes.freezeColumns(1);

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(outputPath);

for (const range of ["A1:V6", `A${Math.max(2, rowCount - 2)}:V${rowCount + 1}`]) {
  const check = await workbook.inspect({
    kind: "region",
    sheetId: "Sheet1",
    range,
    maxChars: 8000,
    tableMaxRows: 8,
    tableMaxCols: 22,
  });
  console.log(check.ndjson);
}
const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(formulaErrors.ndjson);

await fs.mkdir("outputs/previews", { recursive: true });
const preview = await workbook.render({
  sheetName: "Sheet1",
  range: "A1:V26",
  scale: 1.25,
  format: "png",
});
await fs.writeFile(
  "outputs/previews/problem3-Sheet1.png",
  new Uint8Array(await preview.arrayBuffer()),
);

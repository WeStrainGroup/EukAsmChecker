const assert=require("node:assert/strict");
const fs=require("node:fs");
const path=require("node:path");
const crypto=require("node:crypto");
const vm = require("node:vm");
const root=path.resolve(__dirname,".."),out=process.argv[2];
global.sha256=require(path.join(root,"web/dist/vendor/sha256.min.js"));
require(path.join(root,"web/dist/core.js"));
require(path.join(root,"web/dist/data.js"));
require(path.join(root,"web/dist/zip.js"));
// Files are assigned by the chosen input type, with samples matched by name.
const appSource = fs.readFileSync(path.join(root, "web/src/app.js"), "utf8");
const controls = Object.fromEntries(["fastaFiles", "coverageFiles", "ecmFiles"].map(id => [id, {}]));
const inputContext = {
  files: {fasta: [], coverage: [], ecm: []}, samples: [],
  $: id => controls[id], renderSamples() {}, refreshColumns() {}
};
vm.createContext(inputContext);
vm.runInContext(appSource.slice(appSource.indexOf("function path(file)"), appSource.indexOf("function renderSamples()")), inputContext);
controls.fastaFiles.onchange({target: {files: [{name: "a.fasta"}, {name: "b.fasta"}]}});
controls.coverageFiles.onchange({target: {files: [{name: "b.coverage.tsv"}, {name: "a.coverage.tsv"}]}});
controls.ecmFiles.onchange({target: {files: [{name: "a.classification.tsv"}, {name: "b.classification.tsv"}]}});
assert.equal(inputContext.samples.length, 2);
for (const sample of inputContext.samples) {
  assert.equal(inputContext.files.coverage[sample.coverage].name, sample.id + ".coverage.tsv");
  assert.equal(inputContext.files.ecm[sample.ecm].name, sample.id + ".classification.tsv");
}

// Mixed-label caches must preserve filtering, seeded order and threshold-drag reuse.
const reviewSource = fs.readFileSync(path.join(root, "src/eukasmchecker/review.html"), "utf8");
let rasterizedPoints = 0;
const renderRows = [
  {length: 499, coverage: 10, label: "NA"},
  {length: 500, coverage: 10, label: "eukaryota"},
  {length: 1000, coverage: 20, label: "other"},
  {length: 2000, coverage: 30, label: "eukaryota"}
];
function mockCanvas() {
  const canvas = {points: []};
  const context = {
    scale() {}, beginPath() {},
    clearRect() {canvas.points = [];},
    arc(x, y) {this.point = {x, y}; rasterizedPoints++;},
    fill() {canvas.points.push({...this.point, color: this.fillStyle, alpha: this.globalAlpha});}
  };
  canvas.getContext = () => context;
  return canvas;
}
const renderContext = {
  pointLayers: null, hoverPoints: [], rows: renderRows,
  availableLabels: ["eukaryota", "other", "NA"], labelColor: index => String(index),
  document: {createElement: mockCanvas}
};
vm.createContext(renderContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("  function preparePointLayers("), reviewSource.indexOf("  function drawDistributions(")), renderContext);
const geometry = {w: 300, h: 200, left: 20, right: 280, top: 20, bottom: 180,
  x: {a: 2, b: 4, f: length => 20 + (Math.log10(length) - 2) / 2 * 260},
  y: {b: 5, f: coverage => 180 - Math.asinh(coverage) / 5 * 160}};
for (const retain of [true, false]) for (const labels of [["eukaryota"], ["other", "NA"], ["eukaryota", "other", "NA"]]) {
  const settings = {length: 500, coverage: 10, maxCoverage: 20, retain_sample: retain, labels};
  renderContext.result = EACCore.evaluate(renderRows, settings);
  renderContext.preparePointLayers(geometry, 1, settings);
  const cached = renderContext.pointLayers;
  const rendered = rasterizedPoints;
  const shown = renderRows.map(() => []);
  let rectangles = [];
  const context = {
    save() {}, restore() {}, beginPath() {rectangles = [];},
    rect(x, y, w, h) {rectangles.push({x, y, w, h});}, clip() {},
    drawImage(canvas) {
      for (const point of canvas.points) {
        const index = renderContext.hoverPoints.findIndex(item => item.x === point.x && item.y === point.y);
        const inside = rectangles.filter(r => point.x >= r.x && point.x <= r.x + r.w && point.y >= r.y && point.y <= r.y + r.h).length % 2;
        if (inside) shown[index].push(point.alpha);
      }
    }
  };
  renderContext.drawPoints(context, geometry, settings, geometry.x.f(500), geometry.y.f(10), geometry.y.f(20));
  shown.forEach((alphas, index) => assert.deepEqual(alphas, [renderContext.result.annotated[index].keep ? .55 : .22]));
  const changedThresholds = {...settings, length: 800, coverage: 15, maxCoverage: 25};
  renderContext.result = EACCore.evaluate(renderRows, changedThresholds);
  renderContext.preparePointLayers(geometry, 1, changedThresholds);
  assert.equal(renderContext.pointLayers, cached);
  assert.equal(rasterizedPoints, rendered);
}
// Fully overlapping, initially label-sorted input must be interleaved on both layers.
const overlapRows = ["eukaryota", "other", "NA"].flatMap(label =>
  Array.from({length: 30}, () => ({length: 1000, coverage: 20, label})));
renderContext.rows = overlapRows;
renderContext.pointLayers = null;
const allLabels = {length: 500, coverage: 10, maxCoverage: null, retain_sample: true, labels: ["eukaryota", "other", "NA"]};
renderContext.result = EACCore.evaluate(overlapRows, allLabels);
renderContext.preparePointLayers(geometry, 1, allLabels);
const firstOrder = Array.from(renderContext.pointLayers.order);
assert.equal(new Set(firstOrder).size, overlapRows.length);
const colors = renderContext.pointLayers.layers[1].canvas.points.map(point => point.color);
assert.deepEqual(renderContext.pointLayers.layers[0].canvas.points.map(point => point.color), colors);
assert(colors.filter((color, index) => index > 0 && color !== colors[index - 1]).length > 30);
for (const color of ["0", "1", "2"]) assert.equal(colors.filter(value => value === color).length, 30);
const fadedCanvas = renderContext.pointLayers.layers[0].canvas;
const beforeSelection = rasterizedPoints;
const eukOnly = {...allLabels, labels: ["eukaryota"]};
renderContext.result = EACCore.evaluate(overlapRows, eukOnly);
renderContext.preparePointLayers(geometry, 1, eukOnly);
assert.equal(renderContext.pointLayers.layers[0].canvas, fadedCanvas);
assert.equal(rasterizedPoints - beforeSelection, overlapRows.length);
assert.deepEqual(renderContext.pointLayers.layers[1].canvas.points.map(point => point.color), colors);
for (const point of renderContext.pointLayers.layers[1].canvas.points) assert.equal(point.alpha, point.color === "0" ? .55 : .22);
renderContext.preparePointLayers({...geometry, w: 320}, 2, eukOnly);
assert.deepEqual(Array.from(renderContext.pointLayers.order), firstOrder);
renderContext.pointLayers = null;
renderContext.preparePointLayers(geometry, 1, eukOnly);
assert.deepEqual(Array.from(renderContext.pointLayers.order), firstOrder);
renderContext.rows = [];
renderContext.pointLayers = null;
renderContext.result = null;
renderContext.preparePointLayers(geometry, 1, allLabels);
assert.equal(renderContext.pointLayers.layers[0].canvas.points.length, 0);
console.log("Seeded mixed-label order, cache reuse and filtering appearance passed.");

(async()=>{
  const fasta=Buffer.from(">a full header\r\n"+"A".repeat(10001)+"\r\n>b\n"+"C".repeat(500)+"\n>c\n"+"G".repeat(499)+"\n>d tail\r\n"+"T".repeat(935)+"\r\n");
  const coverage="identifier\tmean_reads\nc\t0\nb\t10\nd\t84.9343\na\t100.5\n";
  const ecm="id\tbp\tkingdom\nd\t935\tother\na\t10001\teuk\nc\t499\tNA\nb\t500\tother\n";
  fs.mkdirSync(out,{recursive:true});
  fs.writeFileSync(path.join(out,"assembly.fasta"),fasta);
  fs.writeFileSync(path.join(out,"coverage.tsv"),coverage);
  fs.writeFileSync(path.join(out,"classification.tsv"),ecm);
  const original=new File([fasta],"assembly.fasta");
  const sample={id:"sample_A",fasta:original,coverage:new File([coverage],"coverage.tsv"),ecm:new File([ecm],"classification.tsv")};
  const columns={coverageId:"identifier",coverageValue:"mean_reads",ecmId:"id",ecmLabel:"kingdom",ecmLength:"bp"};
  const prepared=await EACData.prepare(sample,columns);
  const raw=JSON.parse(await prepared.blob.text());
  assert.equal(raw.input_sha256.fasta,crypto.createHash("sha256").update(fasta).digest("hex"));
  const cases=[
    {name:"default",min:10,max:null,labels:["eukaryota","other","NA"],retain:true,count:3,bp:11436},
    {name:"fractional",min:85,max:null,labels:["eukaryota","other","NA"],retain:true,count:1,bp:10001},
    {name:"upper",min:10,max:100,labels:["eukaryota","other","NA"],retain:true,count:2,bp:1435},
    {name:"labels",min:10,max:null,labels:["other"],retain:true,count:2,bp:1435},
    {name:"excluded",min:null,max:null,labels:["eukaryota","other","NA"],retain:false,count:0,bp:0}
  ];
  for(const item of cases){
    const decision={sample_id:"sample_A",fingerprint:prepared.meta.fingerprint,eac_path:prepared.meta.eac_path,
      eac_sha256:prepared.meta.eac_sha256,reviewed:true,retain_sample:item.retain,min_length_bp:500,min_coverage:item.min,
      max_coverage:item.max,labels:item.labels,expected_kept_contigs:item.count,expected_kept_bp:item.bp};
    const result=await EACData.apply(prepared.blob,original,decision,prepared.meta);
    fs.writeFileSync(path.join(out,item.name+".fasta"),Buffer.from(await result.fasta.arrayBuffer()));
    fs.writeFileSync(path.join(out,item.name+".tsv"),await result.table.text());
    const zip=await EACZip.build([{name:"sample_A/sample_A.fasta",blob:result.fasta},{name:"sample_A/sample_A.sequences.tsv",blob:result.table}]);
    fs.writeFileSync(path.join(out,item.name+".zip"),Buffer.from(await zip.arrayBuffer()));
    await assert.rejects(()=>EACData.apply(prepared.blob,original,{...decision,reviewed:false},prepared.meta),/Confirm/);
    await assert.rejects(()=>EACData.apply(prepared.blob,original,{...decision,expected_kept_contigs:999},prepared.meta),/counts/);
  }
  await assert.rejects(()=>EACData.prepare({...sample,coverage:new File([coverage+"b\t10\n"],"duplicate.tsv")},columns),/duplicate/);
  await assert.rejects(()=>EACData.prepare({...sample,coverage:new File([coverage.replace("b\t10\n","")],"missing.tsv")},columns),/same IDs/);
  await assert.rejects(()=>EACData.prepare({...sample,ecm:new File([ecm.replace("b\t500","b\t501")],"length.tsv")},columns),/length differs/);
  await assert.rejects(()=>EACData.prepare({...sample,coverage:new File([coverage.replace("b\t10","b\tNaN")],"nan.tsv")},columns),/Invalid/);
  const noLength=await EACData.prepare(sample,{...columns,ecmLength:""});
  assert.deepEqual(JSON.parse(await noLength.blob.text()).rows,raw.rows);
  const chunks=new File([fasta],"chunked.fasta");
  chunks.stream=()=>new ReadableStream({start(controller){for(let i=0;i<fasta.length;i+=17)controller.enqueue(new Uint8Array(fasta.subarray(i,i+17)));controller.close();}});
  assert.equal((await EACData.fasta(chunks)).hash,raw.input_sha256.fasta);
  fs.writeFileSync(path.join(out,"cases.json"),JSON.stringify(cases));
  console.log("Browser parsing, boundary filters, exclusions, validation and ZIP generation passed.");
})().catch(error=>{console.error(error);process.exitCode=1;});

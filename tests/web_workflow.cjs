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

// Cached point layers must preserve the filter result at every sequence centre.
const reviewSource = fs.readFileSync(path.join(root, "src/eukasmchecker/review.html"), "utf8");
let rasterizedPoints = 0;
const renderRows = [
  {length: 499, coverage: 10, label: "NA"},
  {length: 500, coverage: 10, label: "eukaryota"},
  {length: 1000, coverage: 20, label: "other"},
  {length: 2000, coverage: 30, label: "eukaryota"}
];
const renderContext = {
  pointLayers: null, hoverPoints: [], rows: renderRows,
  availableLabels: ["eukaryota", "other", "NA"], labelColor: () => "black",
  document: {createElement: () => ({getContext: () => ({
    scale() {}, beginPath() {}, arc() {rasterizedPoints++;}, fill() {}
  })})}
};
vm.createContext(renderContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("  function preparePointLayers("), reviewSource.indexOf("  function drawDistributions(")), renderContext);
const geometry = {w: 300, h: 200, left: 20, right: 280, top: 20, bottom: 180,
  x: {a: 2, b: 4, f: length => 20 + (Math.log10(length) - 2) / 2 * 260},
  y: {b: 5, f: coverage => 180 - Math.asinh(coverage) / 5 * 160}};
renderContext.preparePointLayers(geometry, 1);
const originalLayers = renderContext.pointLayers;
renderContext.preparePointLayers(geometry, 1);
assert.equal(renderContext.pointLayers, originalLayers);
assert.equal(rasterizedPoints, renderRows.length * 2);
for (const point of renderContext.hoverPoints) {
  assert.equal(point.xBin, Math.floor((Math.log10(point.r.length) - 2) / 2 * 48));
  assert.equal(point.yBin, Math.floor(Math.asinh(point.r.coverage) / 5 * 48));
}
for (const retain of [true, false]) for (const labels of [["eukaryota"], ["other", "NA"], ["eukaryota", "other", "NA"]]) {
  const state = {length: 500, coverage: 10, maxCoverage: 20, retain_sample: retain, labels};
  renderContext.result = EACCore.evaluate(renderRows, state);
  const shown = renderRows.map(() => []);
  let rectangles = [];
  const context = {
    save() {}, restore() {}, beginPath() {rectangles = [];},
    rect(x, y, w, h) {rectangles.push({x, y, w, h});}, clip() {},
    drawImage(canvas) {
      const label = renderContext.pointLayers.layers.findIndex(pair => pair.some(layer => layer.canvas === canvas));
      const pass = renderContext.pointLayers.layers[label].findIndex(layer => layer.canvas === canvas);
      renderContext.hoverPoints.forEach((point, index) => {
        const inside = rectangles.filter(r => point.x >= r.x && point.x <= r.x + r.w && point.y >= r.y && point.y <= r.y + r.h).length % 2;
        if (point.label === label && inside) shown[index].push(pass);
      });
    }
  };
  renderContext.drawPoints(context, geometry, state, geometry.x.f(500), geometry.y.f(10), geometry.y.f(20));
  shown.forEach((passes, index) => assert.deepEqual(passes, [Number(renderContext.result.annotated[index].keep)]));
}
assert.equal(rasterizedPoints, renderRows.length * 2);
renderContext.preparePointLayers({...geometry, w: 320}, 1);
assert.notEqual(renderContext.pointLayers, originalLayers);
console.log("Point-layer reuse and threshold/label/exclusion rendering passed.");

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

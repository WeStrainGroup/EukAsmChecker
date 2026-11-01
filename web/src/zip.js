/* ZIP64 stored archives preserve FASTA bytes and support archives larger than 4 GiB. */
(function () {
  "use strict";
  const table=new Uint32Array(256);
  for(let n=0;n<256;n++){let c=n;for(let k=0;k<8;k++)c=c&1?0xedb88320^(c>>>1):c>>>1;table[n]=c>>>0;}
  async function crc32(blob) {
    let crc=0xffffffff;
    for await(const chunk of blob.stream())for(const byte of chunk)crc=table[(crc^byte)&255]^(crc>>>8);
    return (crc^0xffffffff)>>>0;
  }
  function pack(fields) {
    const size=fields.reduce((n,[width])=>n+width,0), out=new Uint8Array(size), view=new DataView(out.buffer);
    let at=0;
    for(const [width,value] of fields){
      if(width===8)view.setBigUint64(at,BigInt(value),true);
      else if(width===4)view.setUint32(at,value,true);
      else view.setUint16(at,value,true);
      at+=width;
    }
    return out;
  }
  async function build(entries,progress=()=>{}) {
    const parts=[],central=[];
    let offset=0;
    for(let index=0;index<entries.length;index++){
      const {name,blob}=entries[index];
      if(!name || name.startsWith("/") || name.split("/").includes(".."))throw new Error("Invalid ZIP path.");
      const encoded=new TextEncoder().encode(name);
      if(encoded.length>65535)throw new Error("ZIP filename is too long.");
      const checksum=await crc32(blob),size=blob.size;
      const extra=pack([[2,1],[2,16],[8,size],[8,size]]);
      const header=pack([[4,0x04034b50],[2,45],[2,0x800],[2,0],[2,0],[2,33],[4,checksum],
        [4,0xffffffff],[4,0xffffffff],[2,encoded.length],[2,extra.length]]);
      parts.push(header,encoded,extra,blob);
      const centralExtra=pack([[2,1],[2,24],[8,size],[8,size],[8,offset]]);
      const directory=pack([[4,0x02014b50],[2,45],[2,45],[2,0x800],[2,0],[2,0],[2,33],[4,checksum],
        [4,0xffffffff],[4,0xffffffff],[2,encoded.length],[2,centralExtra.length],[2,0],[2,0],[2,0],[4,0],[4,0xffffffff]]);
      central.push(directory,encoded,centralExtra);
      offset+=header.length+encoded.length+extra.length+size;
      progress(index+1,entries.length);
    }
    const centralSize=central.reduce((n,part)=>n+part.length,0);
    const zip64=pack([[4,0x06064b50],[8,44],[2,45],[2,45],[4,0],[4,0],[8,entries.length],[8,entries.length],[8,centralSize],[8,offset]]);
    const locator=pack([[4,0x07064b50],[4,0],[8,offset+centralSize],[4,1]]);
    const end=pack([[4,0x06054b50],[2,0],[2,0],[2,Math.min(entries.length,65535)],[2,Math.min(entries.length,65535)],
      [4,Math.min(centralSize,0xffffffff)],[4,Math.min(offset,0xffffffff)],[2,0]]);
    return new Blob([...parts,...central,zip64,locator,end],{type:"application/zip"});
  }
  globalThis.EACZip={build,crc32};
})();
